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
        .map_ont()
        .with_index_threads(num_threads)
        .with_sam_out()
        .with_cigar()
        .with_index(ref_path, None)
        .expect("Unable to build index");

    // Set up a channel for passing sequences between threads
    let (sender, receiver): (Sender<(String, Vec<u8>)>, Receiver<(String, Vec<u8>)>) = unbounded();

    // Producer thread: reads sequences and sends them to the channel
    let query_path_clone = query_path.clone();
    let producer = thread::spawn(move || {
        let mut reader = parse_fastx_file(&query_path_clone).expect("valid path/file");
        while let Some(result) = reader.next() {
            let record = result.expect("Error reading record");
            let seq_name = String::from_utf8_lossy(record.id()).into_owned();
            let seq = record.seq().to_vec();

            // Send (sequence_name, sequence_bytes) to the consumers
            sender.send((seq_name, seq)).expect("Error sending data to channel");
        }
    });

    // Consumer threads: receive (name, seq) and perform alignment
    let consumers: Vec<_> = (0..num_threads).map(|_| {
        let receiver = receiver.clone();
        // Make aligner mutable for each thread
        let mut aligner = aligner.clone();
        thread::spawn(move || {
            // For each (name, seq) in the channel, call aligner.map with the name
            let results = receiver
                .iter()
                .filter_map(|(seq_name, seq)| {
                    // Use the query name by passing it as Some(...) in aligner.map
                    aligner.map(&seq, false, false, None, None, Some(seq_name.as_bytes())).ok()
                })
                .collect::<Vec<_>>();
            // Set aligner.idx = None before the thread exits (helps ensure resources free)
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
