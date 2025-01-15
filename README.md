

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

### Pre-built binary (Linux)
```bash
wget https://github.com/jianshu93/adas/releases/download/v0.1.1/adas-Linux-x86-64_v0.1.1.zip
unzip adas-Linux-x86-64_v0.1.0.zip
chmod a+x ./adas-*
./adas-build -h
```
### compile from source
```bash
### Install Rust first if you do not have it
### For Linux
git clone https://github.com/jianshu93/adas
cd adas
cargo build --release --jobs 10
cd target/release


### for MacOS
### Install homebrew first here: https://brew.sh
brew install gcc
git clone https://github.com/jianshu93/adas
cd adas
CC="$(brew --prefix)/bin/gcc-14" CXX="$(brew --prefix)/bin/g++-14" cargo build --release
cd target/release
```

## Usage

1. build HNSW database
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
      --hnsw-ef <HNSW_EF>                     HNSW ef parameter [default: 1600]
      --max_nb_connection <HNSW_MAX_NB_CONN>  HNSW max_nb_conn parameter [default: 256]
      --scale_modify_f <scale_modify>
          scale modification factor in HNSW or HubNSW, must be in [0.2,1] [default: 1.0]
  -h, --help                                  Print help
  -V, --version                               Print version
```
2. Search pre-built HNSW database
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
3. Insert new sequences into HNSW database
```bash
./adas-insert -h
 ************** initializing logger *****************

Insert into Pre-built Hierarchical Navigable Small World Graphs (HNSW) Index

Usage: adas-insert [OPTIONS] --input <FASTA_FILE> --hnsw <DATADIR>

Options:
  -i, --input <FASTA_FILE>  Input FASTA file
  -b, --hnsw <DATADIR>      directory contains pre-built HNSW database files
  -t, --threads <THREADS>   Number of threads for sketching [default: 1]
  -h, --help                Print help
  -V, --version             Print version

```
4. Approximate alignment via chaining
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

## Running test using real-world dataset

```bash
mkdir db_test
cd db_test
adas-build -i ../data/SAR11_cluster_centroid.fa -t 24 -s 256 --hnsw-ef 800 --max_nb_connection 128 --scale_modify_f 0.25
adas-search -i ../data/test_16S_SAR11.fa -b . -n 50 -t 24
```

