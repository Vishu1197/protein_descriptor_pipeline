# Protein Descriptor & Binding Pocket Descriptor Pipeline

A production-grade structural bioinformatics pipeline for extracting comprehensive numerical protein descriptors and binding pocket descriptors from 3D protein structures (.pdb).

Designed for:

- Molecular docking machine learning
- Binding affinity prediction
- Protein-ligand interaction modeling
- Structural bioinformatics
- QSAR/deep learning pipelines
- Virtual screening
- Feature engineering for AI models

---

# Features

## Protein Descriptors

### Sequence-Based Features

- Amino acid composition
- Dipeptide composition
- Molecular weight
- Aromaticity
- Instability index
- Isoelectric point
- GRAVY
- Residue composition statistics

### Structural Features

- Number of atoms
- Number of residues
- Number of chains
- Radius of gyration
- Bounding box dimensions
- Protein volume
- Surface area
- Compactness
- Centroid coordinates

### Secondary Structure Features

- Alpha helix fraction
- Beta sheet fraction
- Coil fraction

### Pocket Descriptors (fpocket)

- Pocket score
- Druggability score
- Pocket volume
- SASA
- Polarity score
- Hydrophobic density
- Alpha sphere statistics

---

# Technologies Used

- Python 3.11
- BioPython
- MDAnalysis
- FreeSASA
- fpocket
- DSSP
- NumPy
- pandas
- tqdm

---

# Installation

---

# Recommended Setup (Windows + WSL)

This repository was developed using:

- Windows 11
- WSL2 Ubuntu
- VS Code Remote WSL
- Conda environment

This setup is highly recommended for structural bioinformatics workflows.

---

# Step 1 — Install WSL

Open PowerShell as Administrator:

```powershell
wsl --install
```

Restart your PC.

Install Ubuntu from Microsoft Store.

Recommended:

- Ubuntu 22.04 LTS

---

# Step 2 — Install Miniconda in WSL

Download:

```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
```

Install:

```bash
bash Miniconda3-latest-Linux-x86_64.sh
```

Initialize conda:

```bash
source ~/.bashrc
```

---

# Step 3 — Create Conda Environment

```bash
conda create -n protein_desc python=3.11 -y
```

Activate:

```bash
conda activate protein_desc
```

---

# Step 4 — Configure Conda Channels

```bash
conda config --add channels conda-forge
conda config --add channels bioconda
conda config --set channel_priority strict
```

---

# Step 5 — Install Dependencies

## Core Scientific Packages

```bash
conda install numpy pandas scipy scikit-learn tqdm biopython -y
```

## MDAnalysis

```bash
conda install mdanalysis -y
```

## FreeSASA

```bash
conda install freesasa -y
```

## DSSP

```bash
conda install dssp -y
```

## fpocket

```bash
conda install fpocket -y
```

## ProDy

```bash
pip install prody
```

---

# Ubuntu Native Installation

If using Ubuntu directly:

```bash
sudo apt update && sudo apt upgrade -y
sudo apt install build-essential git wget curl unzip -y
```

Then follow the same Conda installation steps above.

---

# Verify Installation

## fpocket

```bash
fpocket -h
```

## DSSP

```bash
mkdssp --version
```

## FreeSASA

```bash
python -c "import freesasa; print('FreeSASA OK')"
```

---

# Input Structure

Place all `.pdb` protein structures inside:

```text
proteins/
```

Example:

```text
proteins/
├── EGFR.pdb
├── VEGFR2.pdb
├── 5HT1A.pdb
```

---

# Running the Pipeline

```bash
python protein_descriptor_pipeline.py
```

---

# Output Files

| File | Description |
|------|-------------|
| protein_descriptors.csv | Final ML-ready descriptor dataset |
| failed_proteins.csv | Failed proteins |
| protein_descriptor.log | Runtime logs |

---

# Example Output

| protein | molecular_weight | gravy | radius_of_gyration | pocket_volume |
|---------|------------------|--------|--------------------|---------------|
| EGFR | 134523.2 | -0.32 | 31.2 | 512.4 |

---

# Performance Features

- Multiprocessing support
- Batch processing
- Large-scale protein handling
- Automatic pocket detection
- Logging
- Fault tolerance
- ML-ready numerical outputs

---

# Recommended Hardware

| Component | Recommended |
|-----------|-------------|
| RAM | ≥16 GB |
| CPU | ≥6 cores |
| Storage | SSD |
| OS | Linux / WSL2 |

---

# Scientific Applications

- Docking ML
- Binding affinity prediction
- QSAR
- Protein-ligand interaction modeling
- Structural bioinformatics
- Deep learning feature engineering
- Virtual screening

---

# Citation

If you use this pipeline in research, please cite:

- BioPython
- MDAnalysis
- FreeSASA
- fpocket
- DSSP

---

# License

MIT License

---

# Author

Vishal Chanda

PhD Researcher — Structural Bioinformatics / Computational Biology
