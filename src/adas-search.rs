use clap::{Arg, ArgAction, Command};
use needletail::{parse_fastx_file, Sequence};
use std::path::Path;
use std::sync::{Arc, Mutex};
use std::thread;
use std::path::PathBuf;
use num_cpus;
use std::fs::OpenOptions;

use hnsw_rs::prelude::*;
use gsearch::utils::idsketch::{Id, ItemDict};
use gsearch::utils::parameters::*;
use gsearch::utils::dumpload::*;
use gsearch::utils::reloadhnsw;
use gsearch::utils::SeqDict;
use gsearch::answer::ReqAnswer;

use kmerutils::sketcharg::{SeqSketcherParams, SketchAlgo};
use kmerutils::base::{kmergenerator::*, Kmer32bit, CompressedKmerT};
use kmerutils::sketching::setsketchert::*;
use kmerutils::sketcharg::DataType;
use kmerutils::base::alphabet::Alphabet2b;
use kmerutils::base::sequence::Sequence as SequenceStruct;
use probminhash::setsketcher::SetSketchParams;

use log::{debug, info};
use std::io::BufWriter;
use crossbeam::channel::{unbounded, Sender, Receiver};

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

/// Parameters defining a Request in a Hnsw database
pub struct SearchParams {
    /// directory containing the Hnsw previous dmps
    hnsw_dir : String,
    /// fasta file to search
    search_path : String,
    /// the number of answers by request
    nb_answers : usize,
} // end of RequestParams

impl SearchParams {
    pub fn new(hnsw_dir: String, search_path: String, nb_answers: usize) -> Self {
        SearchParams {
            hnsw_dir,
            search_path,
            nb_answers,
        }
    }

    /// get 
    pub fn get_hnsw_dir(&self) -> &String { &self.hnsw_dir }
    pub fn get_search_path(&self) -> &String { &self.search_path }
    pub fn get_nb_answers(&self) -> usize { self.nb_answers }
}

fn main() {
    // Initialize logger
    println!("\n ************** initializing logger *****************\n");
    let _ = env_logger::Builder::from_default_env().init();

    // Use Clap 4.3 to parse command-line arguments
    let matches = Command::new("adas-search")
        .version("0.1.0")
        .about("Search against Pre-built Hierarchical Navigable Small World Graphs (HNSW) Index")
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
            Arg::new("nb_answers")
                .short('n')
                .long("nbng")
                .value_name("NB_SEARCH_ANSWERS")
                .help("Number of search answers")
                .required(true)
                .action(ArgAction::Set)
                .value_parser(clap::value_parser!(usize))
                .default_value("128"),
        )
        .arg(
            Arg::new("database_path")
                .short('b')
                .long("hnsw")
                .value_name("DATADIR")
                .help("directory contains pre-built HNSW database files")
                .required(true)
                .value_parser(clap::value_parser!(String)),
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
        .get_matches();
    
    let fasta_path = matches.get_one::<String>("input").unwrap().to_string();
    let nb_answers_search = *matches.get_one::<usize>("nb_answers").unwrap();
    let db_path = matches.get_one::<String>("database_path").unwrap().to_string();
    let num_threads = *matches.get_one::<usize>("threads").unwrap();
    println!("Using {} threads", num_threads);
    
    let search_params = SearchParams::new(db_path.clone(), fasta_path.clone(), nb_answers_search);
    let database_dirpath = Path::new(search_params.get_hnsw_dir());

    let hnswio_res = reloadhnsw::get_hnswio(database_dirpath);
    if let Err(e) = hnswio_res {
        panic!("error : {:?}", e);
    }
    let mut hnswio = hnswio_res.unwrap();

    let hnsw_path = std::path::PathBuf::from(database_dirpath);
    let reload_res = ProcessingParams::reload_json(&hnsw_path);
    let processing_params = if let Ok(params) = reload_res {
        log::info!("Sketching parameters: {:?}", params.get_sketching_params());
        log::info!("Block processing: {:?}", params.get_block_flag());
        params
    } else {
        panic!(
            "Cannot reload parameters (file parameters.json) from dir: {:?}",
            &hnsw_path
        );
    };

    
    // Load sketching parameters
    let sketch_params = processing_params.get_sketching_params();

    // Type aliases 
    type Kmer = Kmer32bit;
    type Sketcher = OptDensHashSketch<Kmer, f64>;

    info!("Calling sketch_compressedkmer for OptDensHashSketch::<Kmer32bit, f64>");
    let sketcher = Sketcher::new(&sketch_params);

    println!("Loading HNSW index...");
    // Load the HNSW index
    let hnsw_res = hnswio.load_hnsw::<
        <OptDensHashSketch<Kmer, f64> as SeqSketcherT<Kmer32bit>>::Sig,
        DistHamming
    >();
    if let Err(e) = hnsw_res {
        panic!("error : {:?}", e);
    }
    let hnsw = hnsw_res.unwrap();

    // Load sequence dictionary
    let seqname = database_dirpath.join("seqdict.json");
    log::info!("\n reloading sequence dictionary from {}", &seqname.display());
    let seqdict = SeqDict::reload_json(&seqname);
    let seqdict = match seqdict {
        Ok(seqdict) => seqdict,
        _ => {
            panic!(
                "SeqDict reload from dump file  {} failed",
                seqname.display()
            );
        }
    };
    println!("HNSW index loaded...");

    println!("Sketching...");
    // Set the number of threads globally using Rayon
    rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build_global()
        .unwrap();

    // -------------------------------------------------
    // Use Crossbeam's unbounded channel for multiple consumers
    let (tx, rx): (Sender<(Vec<u8>, Vec<u8>)>, Receiver<(Vec<u8>, Vec<u8>)>) = unbounded();

    // Producer thread: read FASTA -> send to channel
    let fasta_path_clone = fasta_path.clone();
    let producer_handle = thread::spawn(move || {
        let mut reader =
            parse_fastx_file(&Path::new(&fasta_path_clone))
                .expect("Invalid path/file for FASTA");

        while let Some(record) = reader.next() {
            let seqrec = record.expect("Invalid record");
            let seq_id = seqrec.id().to_owned();         // Vec<u8>
            let seq_seq = seqrec.normalize(false).into_owned(); // Vec<u8>
            tx.send((seq_id, seq_seq)).expect("Could not send data");
        }
        // Close sending
        drop(tx);
    });

    // Shared storage
    let signatures = Arc::new(Mutex::new(Vec::new()));
    let itemv = Arc::new(Mutex::new(Vec::new()));

    // We’ll spawn multiple consumer threads
    let mut consumer_handles = Vec::with_capacity(num_threads);

    for _ in 0..num_threads {
        let rx_clone = rx.clone(); // each consumer gets a cloned receiver
        let sketcher_clone = sketcher.clone();
        let signatures_clone = Arc::clone(&signatures);
        let itemv_clone = Arc::clone(&itemv);
        let fasta_path_clone_for_consumer = fasta_path.clone();

        let handle = thread::spawn(move || {
            // Loop over messages
            for (seq_id, seq_seq) in rx_clone.iter() {
                // Convert the sequence
                let seq = ascii_to_seq(&seq_seq).unwrap();
                let vseq = vec![&seq];

                // Sketch the sequence (1 item)
                let signatures_vec = sketcher_clone.sketch_compressedkmer(&vseq, kmer_hash_fn_32bit);
                let signature = signatures_vec.into_iter().next().unwrap();

                // Store the signature
                {
                    let mut s_lock = signatures_clone.lock().unwrap();
                    s_lock.push(signature);
                }
                // Store sequence metadata
                {
                    let mut itemv_lock = itemv_clone.lock().unwrap();
                    let id = Id::new(
                        &fasta_path_clone_for_consumer,
                        &String::from_utf8(seq_id).unwrap(),
                    );
                    let item = ItemDict::new(id, seq.size());
                    itemv_lock.push(item);
                }
            }
        });
        consumer_handles.push(handle);
    }

    // Wait for producer + all consumers
    producer_handle.join().expect("Producer thread panicked");
    for handle in consumer_handles {
        handle.join().expect("Consumer thread panicked");
    }

    // Retrieve the collected signatures and metadata
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
    println!("Sketching done.");
    // Searching
    
    let ef_search = 5000;
    let out_threshold = 1.0;
    let outname = "adas.neighbors.txt";
    let outpath = PathBuf::from(outname);

    let outfile = OpenOptions::new()
        .write(true)
        .create(true)
        .truncate(true)
        .open(&outpath);

    if let Err(e) = outfile {
        log::error!("Could not open file {:?}. Error: {:?}", outpath.as_os_str(), e);
        println!("SeqDict dump: could not open file {:?}", outpath.as_os_str());
        // Early return or panic as needed:
        std::process::exit(1);
    }
    let mut outfile = BufWriter::new(outfile.unwrap());

    println!("Searching HNSW index...");
    // We do parallel_search with our signature vector
    let knn_neighbours = hnsw.parallel_search(&signatures, nb_answers_search, ef_search); 
    for i in 0..knn_neighbours.len() {
        let answer = ReqAnswer::new(i, itemv[i].clone(), &knn_neighbours[i]);
        if answer.dump(&seqdict, out_threshold, &mut outfile).is_err() {
            log::info!(
                "could not dump answer for request id {}",
                answer.get_request_id().get_id().get_fasta_id()
            );
        }
    }
    println!("Searching HNSW index done. Search results saved to adas.neighbors.txt");
}
