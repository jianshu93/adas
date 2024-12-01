use clap::{Arg, ArgAction, Command};
use gsearch::utils::reloadhnsw;
use log::{debug, info};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;
use hnsw_rs::prelude::*;
use annembed_gsearch::fromhnsw::kgraph::KGraph;
use annembed_gsearch::fromhnsw::kgraph_from_hnsw_all;
use num::Float;
use kmerutils::sketcharg::{SeqSketcherParams, SketchAlgo};
use kmerutils::base::{kmergenerator::*, Kmer32bit, CompressedKmerT};
use kmerutils::sketching::setsketchert::*;
use kmerutils::sketcharg::DataType;
use kmerutils::base::alphabet::Alphabet2b;
use kmerutils::base::sequence::Sequence as SequenceStruct;
use num_traits::cast::FromPrimitive;

fn main() {
    // Initialize logger
    println!("\n ************** initializing logger *****************\n");
    let _ = env_logger::Builder::from_default_env().init();

    // Use Clap to parse command-line arguments
    let matches = Command::new("adas-knn")
        .version("0.1.0")
        .about("Extract K Nearest Neighbors (K-NN) from HNSW graph")
        .arg(
            Arg::new("database_path")
                .short('b')
                .long("hnsw")
                .value_name("DATADIR")
                .help("Directory containing pre-built HNSW database files")
                .required(true)
                .value_parser(clap::value_parser!(String)),
        )
        .arg(
            Arg::new("output")
                .short('o')
                .long("output")
                .value_name("OUTPUT_PATH")
                .help("Output path to write the neighbor list")
                .required(true)
                .action(ArgAction::Set)
                .value_parser(clap::value_parser!(String)),
        )
        .get_matches();

    let db_path = matches
        .get_one::<String>("database_path")
        .unwrap()
        .to_string();
    let out_path = matches.get_one::<String>("output").unwrap().to_string();

    let database_dirpath = Path::new(&db_path);
    let hnswio_res = reloadhnsw::get_hnswio(database_dirpath);
    if hnswio_res.is_err() {
        panic!("Error: {:?}", hnswio_res.err());
    }
    let mut hnswio = hnswio_res.unwrap();

    // Load the HNSW graph
    let hnsw_res = hnswio.load_hnsw::<
        <OptDensHashSketch<Kmer32bit, f64> as SeqSketcherT<Kmer32bit>>::Sig,
        DistHamming,
    >();

    if hnsw_res.is_err() {
        panic!("Error: {:?}", hnsw_res.err());
    }
    let hnsw = hnsw_res.unwrap();

    let knbn = 64;

    // Let the compiler infer T and D, specify F as f64
    let kgraph_res = kgraph_from_hnsw_all::<_, _, f64>(&hnsw, knbn);

    match kgraph_res {
        Ok(kgraph) => {
            println!(
                "KGraph successfully created with {} nodes.",
                kgraph.get_nb_nodes()
            );
            // Save the neighbor list to a file
            if let Err(e) = save_neighbor_list_to_file(&kgraph, &out_path) {
                eprintln!("Error saving neighbor list: {:?}", e);
            } else {
                println!("Neighbor list saved to {}", out_path);
            }
        }
        Err(e) => {
            eprintln!("Error creating KGraph: {:?}", e);
        }
    }
}

// Function to save neighbor lists to a file
fn save_neighbor_list_to_file<F>(
    kgraph: &KGraph<F>,
    output_file: &str,
) -> std::io::Result<()>
where
    F: FromPrimitive + Float + std::fmt::UpperExp + Sync + Send + std::iter::Sum,
{
    let file = File::create(output_file)?;
    let mut writer = BufWriter::new(file);

    // Iterate over each node index in the KGraph
    for node_idx in 0..kgraph.get_nb_nodes() {
        // Get the DataId for the current node
        if let Some(node_data_id) = kgraph.get_data_id_from_idx(node_idx) {
            // Get the list of outgoing edges (neighbors)
            let edges = kgraph.get_out_edges_by_idx(node_idx);
            // Write the node's DataId
            write!(writer, "{}:", node_data_id)?;
            // Iterate over each neighbor
            for edge in edges {
                // Get the neighbor's index
                let neighbor_idx = edge.node;
                // Get the neighbor's DataId
                if let Some(neighbor_data_id) = kgraph.get_data_id_from_idx(neighbor_idx) {
                    // Write the neighbor's DataId and the edge weight
                    write!(
                        writer,
                        " {}:{:.6}",
                        neighbor_data_id,
                        edge.weight.to_f64().unwrap()
                    )?;
                }
            }
            // Newline after each node's neighbor list
            writeln!(writer)?;
        }
    }

    writer.flush()?;
    Ok(())
}
