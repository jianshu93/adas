use clap::{Arg, ArgAction, Command};
use needletail::{parse_fastx_file, Sequence};
use std::path::Path;
use std::sync::mpsc::{channel, Sender, Receiver};
use std::sync::{Arc, Mutex};
use std::thread;
use std::path::PathBuf;
use num_cpus;

use hnsw_rs::prelude::*;
use gsearch::utils::idsketch::{Id, ItemDict};
use gsearch::utils::parameters::*;
use gsearch::utils::dumpload::*;
use gsearch::utils::SeqDict;
use kmerutils::sketcharg::{SeqSketcherParams, SketchAlgo};
use kmerutils::base::{kmergenerator::*, Kmer32bit, CompressedKmerT};
use kmerutils::sketching::setsketchert::*;
use kmerutils::sketcharg::DataType;
use kmerutils::base::alphabet::Alphabet2b;
use kmerutils::base::sequence::Sequence as SequenceStruct;
use probminhash::setsketcher::SetSketchParams;
use log::{debug, info};

fn ascii_to_seq(bases: &[u8]) -> Result<SequenceStruct, ()> {
    let alphabet = Alphabet2b::new();
    let mut seq = SequenceStruct::with_capacity(2, bases.len());
    seq.encode_and_add(bases, &alphabet);
    Ok(seq)
} // end of ascii_to_seq

// Define the k-mer hash function as a function
fn kmer_hash_fn_32bit(kmer: &Kmer32bit) -> <Kmer32bit as CompressedKmerT>::Val {
    let canonical = kmer.reverse_complement().min(*kmer);
    let nb_alphabet_bits = Alphabet2b::new().get_nb_bits();
    let mask: <Kmer32bit as CompressedKmerT>::Val =
        ((1u64 << (nb_alphabet_bits * kmer.get_nb_base())) - 1)
            .try_into()
            .unwrap();
    let hashval = canonical.get_compressed_value() & mask;
    hashval
}

fn main() {
    // Initialize logger
    println!("\n ************** initializing logger *****************\n");
    let _ = env_logger::Builder::from_default_env().init();

    // Use Clap 4.3 to parse command-line arguments
    let matches = Command::new("adas-build")
        .version("0.1.0")
        .about("Build Hierarchical Navigable Small World Graphs (HNSW) with MinHash sketching")
        .arg(
            Arg::new("input")
                .short('i')
                .long("input")
                .value_name("FASTA_FILE")
                .help("Input FASTA file")
                .required(true)
                .action(ArgAction::Set)
                .value_parser(clap::value_parser!(String)),
        )
        .arg(
            Arg::new("kmer_size")
                .short('k')
                .long("kmer-size")
                .value_name("KMER_SIZE")
                .help("Size of k-mers, must be ≤14")
                .action(ArgAction::Set)
                .value_parser(clap::value_parser!(usize))
                .default_value("8"),
        )
        .arg(
            Arg::new("sketch_size")
                .short('s')
                .long("sketch-size")
                .value_name("SKETCH_SIZE")
                .help("Size of the sketch")
                .action(ArgAction::Set)
                .value_parser(clap::value_parser!(usize))
                .default_value("512"),
        )
        .arg(
            Arg::new("threads")
                .short('t')
                .long("threads")
                .value_name("THREADS")
                .help("Number of threads for sketching")
                .action(ArgAction::Set)
                .value_parser(clap::value_parser!(usize))
                .default_value("1"),
        )
        .arg(
            Arg::new("hnsw_capacity")
                .long("hnsw-capacity")
                .value_name("HNSW_CAPACITY")
                .help("HNSW capacity parameter")
                .action(ArgAction::Set)
                .value_parser(clap::value_parser!(usize))
                .default_value("50000000"),
        )
        .arg(
            Arg::new("hnsw_ef")
                .long("hnsw-ef")
                .value_name("HNSW_EF")
                .help("HNSW ef parameter")
                .action(ArgAction::Set)
                .value_parser(clap::value_parser!(usize))
                .default_value("1600"),
        )
        .arg(
            Arg::new("hnsw_max_nb_conn")
                .long("max_nb_connection")
                .value_name("HNSW_MAX_NB_CONN")
                .help("HNSW max_nb_conn parameter")
                .action(ArgAction::Set)
                .value_parser(clap::value_parser!(u8))
                .default_value("255"),
        )
        .arg(Arg::new("scale_modification")
            .long("scale_modify_f")
            .help("scale modification factor in HNSW or HubNSW, must be in [0.2,1]")
            .value_name("scale_modify")
            .default_value("1.0")
            .action(ArgAction::Set)
            .value_parser(clap::value_parser!(f64))
        )
        .get_matches();

    // Parse the command-line arguments using get_one()
    let fasta_path = matches.get_one::<String>("input").unwrap().to_string();
    let kmer_size = *matches.get_one::<usize>("kmer_size").unwrap();
    let sketch_size = *matches.get_one::<usize>("sketch_size").unwrap();
    let num_threads = *matches.get_one::<usize>("threads").unwrap();
    println!("Using {} threads", num_threads);
    let hnsw_capacity = *matches.get_one::<usize>("hnsw_capacity").unwrap();
    let hnsw_ef = *matches.get_one::<usize>("hnsw_ef").unwrap();
    let hnsw_max_nb_conn = *matches.get_one::<u8>("hnsw_max_nb_conn").unwrap();
    let scale_modify = *matches.get_one::<f64>("scale_modification").unwrap();
    if kmer_size > 15 {
        panic!("kmer_size must be ≤15");
    }

    // Set the number of threads globally using Rayon
    rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build_global()
        .unwrap();

    // Set up sketching parameters
    let sketch_args = SeqSketcherParams::new(
        kmer_size,
        sketch_size,
        SketchAlgo::OPTDENS,
        DataType::DNA,
    );

    // Define type aliases for clarity
    type Kmer = Kmer32bit;
    type Sketcher = OptDensHashSketch<Kmer, f64>;

    info!("Calling sketch_compressedkmer for OptDensHashSketch::<Kmer32bit, f64>");
    let sketcher = Sketcher::new(&sketch_args);

    // Create a channel for the producer-consumer model
    let (tx, rx): (Sender<(Vec<u8>, Vec<u8>)>, Receiver<(Vec<u8>, Vec<u8>)>) = channel();

    // Spawn a thread to read the FASTA file and send sequences to the channel
    let fasta_path_clone = fasta_path.clone();
    let producer_handle = thread::spawn(move || {
        // Open the FASTA file
        let mut reader =
            parse_fastx_file(&Path::new(&fasta_path_clone)).expect("Invalid path/file");

        // Read sequences and send them to the channel
        while let Some(record) = reader.next() {
            let seqrec = record.expect("Invalid record");
            let seq_id = seqrec.id().to_owned();
            let seq_seq = seqrec.normalize(false).into_owned();
            tx.send((seq_id, seq_seq)).expect("Could not send data");
        }
        // Close the channel when done
        drop(tx);
    });

    // Set up storage for signatures and sequence metadata
    let signatures = Arc::new(Mutex::new(Vec::new()));
    let signatures_clone = Arc::clone(&signatures);
    let itemv = Arc::new(Mutex::new(Vec::new())); // For seqdict
    let itemv_clone = Arc::clone(&itemv);

    // Clone the sketcher for use in the consumer thread
    let sketcher_clone = sketcher.clone();

    let fasta_path_clone_for_consumer = fasta_path.clone();

    // Spawn a thread to receive sequences and sketch them
    let consumer_handle = thread::spawn(move || {
        // Receive sequences and sketch them
        for (seq_id, seq_seq) in rx {
            // Convert the sequence to the format needed by the sketcher
            let seq = ascii_to_seq(&seq_seq).unwrap();
            let vseq = vec![&seq];

            // Sketch the sequence
            let signatures_vec =
                sketcher_clone.sketch_compressedkmer(&vseq, kmer_hash_fn_32bit);

            // Since we have one sequence, signatures_vec has one element
            let signature = signatures_vec.into_iter().next().unwrap();

            // Store the signature
            {
                let mut signatures_lock = signatures_clone.lock().unwrap();
                signatures_lock.push(signature);
            }
            // Store sequence metadata
            {
                let mut itemv_lock = itemv_clone.lock().unwrap();
                // Create an Id instance using the constructor method
                let id = Id::new(
                    &fasta_path_clone_for_consumer.clone(),
                    &String::from_utf8(seq_id.clone()).unwrap(),
                );
                // Create an ItemDict instance
                let item = ItemDict::new(id, seq.size());
                itemv_lock.push(item);
            }
        }
    });

    // Wait for both threads to finish
    producer_handle.join().expect("Producer thread panicked");
    consumer_handle.join().expect("Consumer thread panicked");

    // After collecting all signatures, build the HNSW index in the main thread
    println!("Building HNSW index...");

    // Retrieve the collected signatures and sequence metadata
    let signatures = Arc::try_unwrap(signatures).unwrap().into_inner().unwrap();
    let mut itemv = match Arc::try_unwrap(itemv) {
        Ok(mutex) => mutex.into_inner().unwrap(),
        Err(_arc_mutex) => {
            // Handle the case where the Arc has more than one strong reference.
            panic!("Cannot unwrap Arc because there are multiple references");
        }
    };

    // Create data as Vec<(&Vec<f64>, usize)> for HNSW insertion
    let data: Vec<(&Vec<f64>, usize)> = signatures.iter().enumerate().map(|(idx, sig)| (sig, idx)).collect();

    // Build the HNSW index
    let hnsw_params = HnswParams::new(hnsw_capacity, hnsw_ef, hnsw_max_nb_conn, scale_modify);

    let mut hnsw = Hnsw::<
        <Sketcher as SeqSketcherT<Kmer>>::Sig,
        DistHamming,
    >::new(
        hnsw_params.get_max_nb_connection() as usize,
        hnsw_params.capacity,
        16, // Adjust as needed
        hnsw_params.get_ef(),
        DistHamming {},
    );
    hnsw.modify_level_scale(scale_modify);
    hnsw.set_extend_candidates(true);
    hnsw.set_keeping_pruned(false);
    
    // Parallel insert all signatures to build HNSW index, using all threads by default
    hnsw.parallel_insert(&data);

    // Create seqdict
    let mut seqdict = SeqDict::new(1000000);

    seqdict.0.append(&mut itemv);

    // Now, create processing_params, which includes HnswParams, sketch_args, and block_flag
    let block_flag = false; // Set as appropriate
    let processing_params = ProcessingParams::new(hnsw_params, sketch_args, block_flag);

    // Now dump all data
    let dump_path = PathBuf::from(".");
    let dump_path_ref = &dump_path;

    let _ = dumpall(dump_path_ref, &hnsw, &seqdict, &processing_params);

    println!("HNSW index built successfully \n");

}
