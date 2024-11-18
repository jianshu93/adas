

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

## Usage
```bash
adas-build -h

adas-search -h

```

