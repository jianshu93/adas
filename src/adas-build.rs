use clap::{Arg, ArgAction, Command};
use crossbeam::channel::{unbounded, Sender, Receiver};
use needletail::{parse_fastx_file, Sequence};
use std::path::Path;
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
    // Initialize logger (optional)
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
        .arg(
            Arg::new("scale_modification")
                .long("scale_modify_f")
                .help("scale modification factor in HNSW or HubNSW, must be in [0.2,1]")
                .value_name("scale_modify")
                .default_value("1.0")
                .action(ArgAction::Set)
                .value_parser(clap::value_parser!(f64))
        )
        .get_matches();

    // Parse the command-line arguments
    let fasta_path = matches.get_one::<String>("input").unwrap().to_string();
    let kmer_size = *matches.get_one::<usize>("kmer_size").unwrap();
    let sketch_size = *matches.get_one::<usize>("sketch_size").unwrap();
    let num_threads = *matches.get_one::<usize>("threads").unwrap();
    println!("Using {} threads", num_threads);
    let hnsw_ef = *matches.get_one::<usize>("hnsw_ef").unwrap();
    let hnsw_max_nb_conn = *matches.get_one::<u8>("hnsw_max_nb_conn").unwrap();
    let scale_modify = *matches.get_one::<f64>("scale_modification").unwrap();

    if kmer_size > 15 {
        panic!("kmer_size must be ≤14");
    }

    
    // If your code uses Rayon for something, set up the Rayon thread pool
    rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build_global()
        .unwrap();

    println!("Sketching...");
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

    // Use Crossbeam's unbounded channel:
    let (tx, rx): (Sender<(Vec<u8>, Vec<u8>)>, Receiver<(Vec<u8>, Vec<u8>)>) = unbounded();

    // Spawn a producer thread to read the FASTA file and send sequences
    let fasta_path_clone = fasta_path.clone();
    let producer_handle = thread::spawn(move || {
        let mut reader = parse_fastx_file(Path::new(&fasta_path_clone))
            .expect("Invalid path/file for FASTA");

        // Read sequences and send them to the channel
        while let Some(record) = reader.next() {
            let seqrec = record.expect("Invalid record");
            let seq_id = seqrec.id().to_owned(); // Vec<u8>
            let seq_seq = seqrec.normalize(false).into_owned(); // Vec<u8>
            tx.send((seq_id, seq_seq)).expect("Could not send data");
        }
        // Close the sending side
        drop(tx);
    });

    // Set up shared storage for signatures and metadata
    let signatures = Arc::new(Mutex::new(Vec::new()));
    let itemv = Arc::new(Mutex::new(Vec::new()));

    // We’ll spawn multiple consumer threads
    let mut consumer_handles = Vec::with_capacity(num_threads);

    for _ in 0..num_threads {
        // Clone everything needed inside this thread:
        let rx_clone = rx.clone(); // Crossbeam receivers can be cloned
        let sketcher_clone = sketcher.clone();
        let signatures_clone = Arc::clone(&signatures);
        let itemv_clone = Arc::clone(&itemv);
        let fasta_path_clone_for_consumer = fasta_path.clone();

        let handle = thread::spawn(move || {
            // Each consumer thread pulls data in parallel
            for (seq_id, seq_seq) in rx_clone.iter() {
                // Convert the sequence to the format needed by the sketcher
                let seq = ascii_to_seq(&seq_seq).unwrap();
                let vseq = vec![&seq];

                // Sketch the sequence
                let signatures_vec =
                    sketcher_clone.sketch_compressedkmer(&vseq, kmer_hash_fn_32bit);

                // We only have one sequence in the vector
                let signature = signatures_vec.into_iter().next().unwrap();

                // Store the signature (only 1 item)
                {
                    let mut signatures_lock = signatures_clone.lock().unwrap();
                    signatures_lock.push(signature);
                }
                // Store sequence metadata
                {
                    let mut itemv_lock = itemv_clone.lock().unwrap();
                    // Create an Id instance
                    let id = Id::new(
                        &fasta_path_clone_for_consumer,
                        &String::from_utf8(seq_id).unwrap(),
                    );
                    // Create an ItemDict instance
                    let item = ItemDict::new(id, seq.size());
                    itemv_lock.push(item);
                }
                // Here seq and seq_seq are dropped -> memory freed
            }
        });
        consumer_handles.push(handle);
    }

    // Wait for the producer to finish
    producer_handle.join().expect("Producer thread panicked");

    // Wait for all consumers to finish
    for handle in consumer_handles {
        handle.join().expect("Consumer thread panicked");
    }
    println!("Sketching done...");
    // ---- Now build HNSW in the main thread ----
    println!("Building HNSW index...");

    // Retrieve the collected signatures and sequence metadata
    let signatures = Arc::try_unwrap(signatures)
        .expect("Multiple references to signatures remain")
        .into_inner()
        .expect("Mutex was poisoned in signatures");
    let mut itemv = match Arc::try_unwrap(itemv) {
            Ok(mutex) => mutex.into_inner().unwrap(),
            Err(_arc_mutex) => {
                // Handle the case where the Arc has more than one strong reference.
                panic!("Cannot unwrap Arc because there are multiple references");
            }
        };

    // Create data as Vec<(&Vec<f64>, usize)> for HNSW insertion
    let data: Vec<(&Vec<f64>, usize)> = signatures
        .iter()
        .enumerate()
        .map(|(idx, sig)| (sig, idx))
        .collect();

    // Build the HNSW index
    let max_nb_conn: u8 = 255.min(hnsw_max_nb_conn as u8);
    let hnsw_params = HnswParams::new(2_500_000, hnsw_ef, max_nb_conn, scale_modify);

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

    // Parallel insert all signatures to build HNSW index
    hnsw.parallel_insert(&data);

    // Create seqdict
    let mut seqdict = SeqDict::new(1000000);
    seqdict.0.append(&mut itemv);

    // Create processing_params
    let block_flag = false; // or true, as needed
    let processing_params = ProcessingParams::new(hnsw_params, sketch_args, block_flag);

    // Dump all data
    let dump_path = PathBuf::from(".");
    let dump_path_ref = &dump_path;
    let _ = dumpall(dump_path_ref, &hnsw, &seqdict, &processing_params);

    println!("HNSW index built successfully.\n");
}
