# 🧬 Protein Descriptor & Binding Pocket Descriptor Pipeline

<div align="center">

![Python](https://img.shields.io/badge/Python-3.11-blue?style=for-the-badge&logo=python)
![Bioinformatics](https://img.shields.io/badge/Bioinformatics-Structural%20Biology-green?style=for-the-badge)
![Machine Learning](https://img.shields.io/badge/ML-Docking%20Prediction-orange?style=for-the-badge)
![Platform](https://img.shields.io/badge/Platform-WSL%20%7C%20Linux-red?style=for-the-badge)

🚀 **Production-Grade Structural Bioinformatics Pipeline**  
for extracting **protein descriptors** and **binding pocket descriptors** from 3D protein structures (`.pdb`)

</div>

---

# ✨ Overview

This repository provides a high-performance computational biology pipeline for generating **ML-ready numerical descriptors** from protein structures.

Designed for:

✅ Molecular Docking ML  
✅ Binding Affinity Prediction  
✅ Protein-Ligand Interaction Modeling  
✅ Structural Bioinformatics  
✅ QSAR & Deep Learning  
✅ Virtual Screening  
✅ AI-driven Drug Discovery  

---

# 🔥 Features

# 🧪 Protein Descriptor Extraction

## 📌 Sequence-Based Features

- Amino Acid Composition
- Dipeptide Composition
- Molecular Weight
- Aromaticity
- Instability Index
- Isoelectric Point (pI)
- GRAVY Hydrophobicity
- Residue Statistics

---

## 🏗 Structural Features

- Number of Atoms
- Number of Residues
- Number of Chains
- Radius of Gyration
- Bounding Box Dimensions
- Protein Volume
- Solvent Accessible Surface Area (SASA)
- Compactness
- Centroid Coordinates

---

## 🧬 Secondary Structure Features

- Alpha Helix Fraction
- Beta Sheet Fraction
- Coil Fraction

---

# 🕳 Binding Pocket Descriptors (fpocket)

- Pocket Score
- Druggability Score
- Pocket Volume
- Pocket SASA
- Polarity Score
- Hydrophobic Density
- Alpha Sphere Statistics
- Pocket Compactness

---

# ⚡ Performance Features

✅ Multiprocessing Support  
✅ Batch Protein Processing  
✅ Fault-Tolerant Execution  
✅ Large Dataset Handling  
✅ Automatic Pocket Detection  
✅ Logging & Error Tracking  
✅ ML-Ready Numerical Outputs  

---

# 🛠 Technologies Used

| Tool | Purpose |
|------|----------|
| 🐍 Python 3.11 | Core Programming |
| 🧬 BioPython | Protein Parsing |
| 📊 MDAnalysis | Structural Analysis |
| 🌊 FreeSASA | Surface Area Calculation |
| 🕳 fpocket | Pocket Detection |
| 🧱 DSSP | Secondary Structure |
| 🔢 NumPy | Numerical Computing |
| 🐼 pandas | Data Handling |
| ⚡ tqdm | Progress Bars |

---

# 💻 Recommended Setup

<div align="center">

## 🪟 Windows + WSL2 + VS Code

</div>

This repository was developed using:

- 🪟 Windows 11
- 🐧 WSL2 Ubuntu
- 💻 VS Code Remote WSL
- 🐍 Conda Environment

✅ Highly Recommended for Computational Biology Workflows

---

# 🚀 Installation Guide

# 1️⃣ Install WSL2

Open **PowerShell as Administrator**:

```powershell
wsl --install
```

Restart your PC.

Then install:

- 🐧 Ubuntu 22.04 LTS

from Microsoft Store.

---

# 2️⃣ Install Miniconda Inside WSL

Download installer:

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

# 3️⃣ Create Conda Environment

```bash
conda create -n protein_desc python=3.11 -y
```

Activate environment:

```bash
conda activate protein_desc
```

---

# 4️⃣ Configure Conda Channels

```bash
conda config --add channels conda-forge
conda config --add channels bioconda
conda config --set channel_priority strict
```

---

# 5️⃣ Install Dependencies

## 🔬 Core Scientific Packages

```bash
conda install numpy pandas scipy scikit-learn tqdm biopython -y
```

---

## 📊 MDAnalysis

```bash
conda install mdanalysis -y
```

---

## 🌊 FreeSASA

```bash
conda install freesasa -y
```

---

## 🧱 DSSP

```bash
conda install dssp -y
```

---

## 🕳 fpocket

```bash
conda install fpocket -y
```

---

## 🧬 ProDy

```bash
pip install prody
```

---

# 🐧 Ubuntu Native Installation

If using native Ubuntu:

```bash
sudo apt update && sudo apt upgrade -y
sudo apt install build-essential git wget curl unzip -y
```

Then follow the same Conda setup instructions above.

---

# ✅ Verify Installation

## 🕳 fpocket

```bash
fpocket -h
```

---

## 🧱 DSSP

```bash
mkdssp --version
```

---

## 🌊 FreeSASA

```bash
python -c "import freesasa; print('FreeSASA OK')"
```

---

# 📂 Input Structure

Place all `.pdb` files inside:

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

# ▶ Running the Pipeline

```bash
python protein_descriptor_pipeline.py
```

---

# 📤 Output Files

| 📄 File | 📌 Description |
|---|---|
| `protein_descriptors.csv` | Final ML-ready descriptor dataset |
| `failed_proteins.csv` | Failed proteins |
| `protein_descriptor.log` | Runtime & error logs |

---

# 📊 Example Output

| protein | molecular_weight | gravy | radius_of_gyration | pocket_volume |
|---|---|---|---|---|
| EGFR | 134523.2 | -0.32 | 31.2 | 512.4 |

---

# 🧠 Scientific Applications

✅ Docking ML  
✅ Binding Affinity Prediction  
✅ Protein-Ligand Interaction Modeling  
✅ Structural Bioinformatics  
✅ Deep Learning Feature Engineering  
✅ QSAR Modeling  
✅ Virtual Screening  
✅ AI Drug Discovery  

---

# ⚙ Recommended Hardware

| Component | Recommended |
|---|---|
| 🧠 RAM | ≥16 GB |
| ⚡ CPU | ≥6 Cores |
| 💾 Storage | SSD |
| 🐧 OS | Linux / WSL2 |

---

# 📚 Citation

If you use this repository in research, please cite:

- BioPython
- MDAnalysis
- FreeSASA
- fpocket
- DSSP

---

# 🤝 Contributing

Contributions are welcome!

Feel free to:

- Fork the repository
- Open issues
- Submit pull requests
- Improve descriptor extraction
- Add new structural features

---

# 📜 License

MIT License © 2026

---

# 👨‍🔬 Author

## Vishal Chanda

🧬 PhD Researcher  
🔬 Structural Bioinformatics & Computational Biology  
🚀 AI for Drug Discovery

---

<div align="center">

# ⭐ If you found this project useful, consider starring the repository!

</div>
