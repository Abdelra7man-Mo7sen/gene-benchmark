# Benchmark of Automated Gene Identification Tools

![Python](https://img.shields.io/badge/Python-3.10%2B-blue)
![Nextflow](https://img.shields.io/badge/Nextflow-23.10.0-green)
![License](https://img.shields.io/badge/License-MIT-yellow)

A large-scale empirical benchmark comparing 11 gene prediction tools on 113 genomes using multi-omics validation.

## Key Results
- **Prokaryotic (complete)**: Prodigal F1=0.987
- **Fragmented metagenomes**: MetaGeneAnnotator F1=0.914  
- **Eukaryotic with RNA-seq**: BRAKER2 F1=0.901
- **Ensemble (XGBoost)**: F1=0.996 / 0.934

## Quick Start
```bash
# Clone repository
git clone https://github.com/Abdelra7man-Mo7sen/gene-benchmark.git
cd gene-benchmark

# Run with Docker
docker pull gene_benchmark/v1.0
docker run -v $(pwd)/data:/data gene_benchmark/v1.0

# Or with Nextflow
nextflow run main.nf -profile test


## Repository Structure
```
├── data/           # Genome datasets (download via script)
├── tools/          # 11 gene prediction tools (Dockerized)
├── scripts/        # Python/R analysis scripts
├── results/        # Output tables and figures
└── docs/           # Supplementary materials
```

## Results Summary Table
| Tool | Prokaryote F1 | Eukaryote F1 | Time (min/Mbp) |
|------|---------------|---------------|----------------|
| Prodigal | 0.987 | — | 0.04 |
| BRAKER2 | — | 0.901 | 187 |
| DeepGene | 0.957 | 0.744 | 0.67 |
| MetaGeneAnnotator | 0.904 | — | 0.11 |

## Citation
If you use this benchmark, please cite:
```
[Author], et al. (2026). Large-scale empirical benchmark of automated gene identification tools. Bioinformatics, 42(5), 1234-1245.
```

## License
MIT - see [LICENSE](LICENSE) file