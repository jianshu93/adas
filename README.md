

# ADAS:Advanced Database Search for Long Sequences
This crate (currently in development) is designed for searching long sequence databases, especially biological sequences. The core idea is from innovative applications of MinHash, Coreset and hierarchical navigable small world graphs (HNSW) algorithms.

<div align="center">
  <img width="75%" src ="adas.jpg">
</div>

## Schematic overview of ADAS
Here we simply describe the algorithm:
1. build HNSW graph for all sequences as a database, the underlying algorithm will be Order MinHash, an LSH for Edit distance.
2. Search new sequences against the pre-built database. The same Order MinHash algorithm will be used.
3. For closest sequences to query seqeunces (e.g., top 100), seed-chain-exend and adaptive banded dynamic progrmming will be performed. 

## Install
```bash
### Install Rust first if you do not have it
git clone https://github.com/jianshu93/adas
cargo build --release

```

## Usage
```bash
adas-build -h

 ************** initializing logger *****************

MinHash sketching and Hierarchical Navigable Small World Graphs (HNSW) building for Long Sequences

Usage: adas-build [OPTIONS] --input <FASTA_FILE>

Options:
  -i, --input <FASTA_FILE>                    Input FASTA file
  -k, --kmer-size <KMER_SIZE>                 Size of k-mers, must be ≤14 [default: 8]
  -s, --sketch-size <SKETCH_SIZE>             Size of the sketch [default: 512]
  -t, --threads <THREADS>                     Number of threads for sketching [default: 1]
      --hnsw-capacity <HNSW_CAPACITY>         HNSW capacity parameter [default: 50000000]
      --hnsw-ef <HNSW_EF>                     HNSW ef parameter [default: 1600]
      --max_nb_connection <HNSW_MAX_NB_CONN>  HNSW max_nb_conn parameter [default: 256]
  -h, --help                                  Print help
  -V, --version                               Print version
```

```bash
adas-search -h

 ************** initializing logger *****************

Search Query Sequences against Pre-built Hierarchical Navigable Small World Graphs (HNSW) Index

Usage: adas-search [OPTIONS] --input <FASTA_FILE> --nbng <NB_SEARCH_ANSWERS> --hnsw <DATADIR>

Options:
  -i, --input <FASTA_FILE>        Input FASTA file
  -n, --nbng <NB_SEARCH_ANSWERS>  Number of search answers [default: 128]
  -b, --hnsw <DATADIR>            directory contains pre-built HNSW database files
  -t, --threads <THREADS>         Number of threads for sketching [default: 1]
  -h, --help                      Print help
  -V, --version                   Print version
```
```bash
adas-chain -h
Long Reads Alignment via Anchor Chaining

Usage: adas-chain [OPTIONS] --reference <REFERENCE_FASTA> --query <QUERY_FASTA> --output <OUTPUT_PATH>

Options:
  -r, --reference <REFERENCE_FASTA>  Reference FASTA file
  -q, --query <QUERY_FASTA>          Query FASTA file
  -t, --threads <THREADS>            Number of threads (default 1) [default: 1]
  -o, --output <OUTPUT_PATH>         Output path to write the results
  -h, --help                         Print help
  -V, --version                      Print version
```

