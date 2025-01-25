use minimap2::Aligner;
use needletail::parse_fastx_file;
use crossbeam_channel::{unbounded, Receiver, Sender};
use std::thread;
use clap::{Arg, ArgAction, Command, value_parser};
use std::fs::File;
use std::io::Write;
use num_cpus;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // Parse command-line arguments using Clap 4.3
    let matches = Command::new("adas-chaining")
        .version("0.1.1")
        .about("Long Reads Alignment via Anchor Chaining")
        .arg(
            Arg::new("reference")
                .short('r')
                .long("reference")
                .value_name("REFERENCE_FASTA")
                .help("Reference FASTA file")
                .required(true)
                .action(ArgAction::Set)
                .value_parser(value_parser!(String)),
        )
        .arg(
            Arg::new("query")
                .short('q')
                .long("query")
                .value_name("QUERY_FASTA")
                .help("Query FASTA file")
                .required(true)
                .action(ArgAction::Set)
                .value_parser(value_parser!(String)),
        )
        .arg(
            Arg::new("threads")
                .short('t')
                .long("threads")
                .value_name("THREADS")
                .help("Number of threads (default 1)")
                .action(ArgAction::Set)
                .value_parser(value_parser!(usize))
                .default_value("1"),
        )
        .arg(
            Arg::new("output")
                .short('o')
                .long("output")
                .value_name("OUTPUT_PATH")
                .help("Output path to write the results (PAF format)")
                .required(true)
                .action(ArgAction::Set)
                .value_parser(value_parser!(String)),
        )
        .get_matches();

    let ref_path = matches.get_one::<String>("reference").unwrap();
    let query_path = matches.get_one::<String>("query").unwrap();
    let num_threads = *matches.get_one::<usize>("threads").unwrap();
    let output_path = matches.get_one::<String>("output").unwrap();

    let num_cpus = num_cpus::get();
    let num_threads = if num_threads > num_cpus {
        num_cpus
    } else {
        num_threads
    };
    // Build the aligner with the specified reference and number of threads
    let aligner = Aligner::builder()
        .ava_ont()
        .with_index_threads(num_threads)
        .with_index(ref_path, None)
        .expect("Unable to build index");

    // Set up a channel for passing sequences between threads
    let (sender, receiver): (Sender<Vec<u8>>, Receiver<Vec<u8>>) = unbounded();

    // Producer thread: reads sequences and sends them to the channel
    let query_path_clone = query_path.clone();
    let producer = thread::spawn(move || {
        let mut reader = parse_fastx_file(&query_path_clone).expect("valid path/file");
        while let Some(result) = reader.next() {
            let record = result.expect("Error reading record");
            sender.send(record.seq().to_vec()).expect("Error sending sequence");
        }
    });

    // Consumer threads: receive sequences and perform alignment
    let consumers: Vec<_> = (0..num_threads).map(|_| {
        let receiver = receiver.clone();
        // Make aligner mutable
        let mut aligner = aligner.clone(); 
        thread::spawn(move || {
            let results = receiver.iter().filter_map(|seq: Vec<u8>| {
                aligner.map(&seq, false, false, None, None, None).ok()
            }).collect::<Vec<_>>();
    
            // Set aligner.idx to None before the thread exits
            aligner.idx = None;
            results
        })
    }).collect();

    // Wait for the producer to finish reading
    producer.join().expect("Producer thread panicked");

    // Collect results from consumer threads
    let mut results = Vec::new();
    for consumer in consumers {
        let mut res = consumer.join().expect("Consumer thread panicked");
        results.append(&mut res);
    }

    // Write the results to output file
    let mut output_file = File::create(output_path)?;
    for result in results {
        writeln!(output_file, "{:?}", result)?;
    }

    Ok(())
}
