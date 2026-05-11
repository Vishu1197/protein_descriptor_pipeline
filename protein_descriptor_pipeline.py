import os
import re
import sys
import math
import glob
import time
import shutil
import logging
import warnings
import traceback
import subprocess
from pathlib import Path
from multiprocessing import Pool, cpu_count

import numpy as np
import pandas as pd
from tqdm import tqdm

from Bio.PDB import PDBParser, DSSP, PPBuilder
from Bio.SeqUtils.ProtParam import ProteinAnalysis

import freesasa
import MDAnalysis as mda
from MDAnalysis.analysis import rms

warnings.filterwarnings("ignore")

# =====================================================================================
# CONFIGURATION
# =====================================================================================

INPUT_DIR = "/home/vishu/dock_predict_model/proteins"
OUTPUT_CSV = "protein_descriptors.csv"
FAILED_CSV = "failed_proteins.csv"
LOG_FILE = "protein_descriptor.log"

USE_MULTIPROCESSING = True
N_CORES = max(1, cpu_count() - 1)

# fpocket executable
# IMPORTANT:
# Add fpocket executable path here if not in PATH
FPOCKET_EXECUTABLE = "fpocket"

# =====================================================================================
# LOGGING
# =====================================================================================

logging.basicConfig(
    filename=LOG_FILE,
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s"
)

console = logging.StreamHandler()
console.setLevel(logging.INFO)
logging.getLogger("").addHandler(console)

# =====================================================================================
# AMINO ACID GROUPS
# =====================================================================================

HYDROPHOBIC = set(["A", "V", "I", "L", "M", "F", "W", "Y"])
POLAR = set(["S", "T", "N", "Q", "C"])
POSITIVE = set(["K", "R", "H"])
NEGATIVE = set(["D", "E"])
AROMATIC = set(["F", "W", "Y", "H"])

# =====================================================================================
# UTILITY FUNCTIONS
# =====================================================================================

def safe_div(a, b):
    if b == 0:
        return 0
    return a / b


def compute_radius_of_gyration(coords):
    centroid = np.mean(coords, axis=0)
    rg = np.sqrt(np.mean(np.sum((coords - centroid) ** 2, axis=1)))
    return rg


def compute_bounding_box(coords):
    mins = coords.min(axis=0)
    maxs = coords.max(axis=0)
    dims = maxs - mins
    return dims


def compute_compactness(volume, surface_area):
    try:
        return (surface_area ** 3) / (36 * np.pi * (volume ** 2))
    except:
        return np.nan


def get_sequence_from_structure(structure):
    ppb = PPBuilder()
    seq = ""
    for pp in ppb.build_peptides(structure):
        seq += str(pp.get_sequence())
    return seq


def residue_statistics(sequence):
    total = len(sequence)

    hydrophobic = sum(aa in HYDROPHOBIC for aa in sequence)
    polar = sum(aa in POLAR for aa in sequence)
    positive = sum(aa in POSITIVE for aa in sequence)
    negative = sum(aa in NEGATIVE for aa in sequence)
    aromatic = sum(aa in AROMATIC for aa in sequence)

    return {
        "hydrophobic_residue_fraction": safe_div(hydrophobic, total),
        "polar_residue_fraction": safe_div(polar, total),
        "positive_residue_fraction": safe_div(positive, total),
        "negative_residue_fraction": safe_div(negative, total),
        "aromatic_residue_fraction": safe_div(aromatic, total),
    }


# =====================================================================================
# SEQUENCE DESCRIPTORS
# =====================================================================================

def calculate_sequence_descriptors(sequence):

    descriptors = {}

    if len(sequence) == 0:
        return descriptors

    analysis = ProteinAnalysis(sequence)

    descriptors["molecular_weight"] = analysis.molecular_weight()
    descriptors["aromaticity"] = analysis.aromaticity()
    descriptors["instability_index"] = analysis.instability_index()
    descriptors["isoelectric_point"] = analysis.isoelectric_point()
    descriptors["gravy"] = analysis.gravy()

    try:
        helix, turn, sheet = analysis.secondary_structure_fraction()
        descriptors["helix_fraction_seq"] = helix
        descriptors["turn_fraction_seq"] = turn
        descriptors["sheet_fraction_seq"] = sheet
    except:
        descriptors["helix_fraction_seq"] = np.nan
        descriptors["turn_fraction_seq"] = np.nan
        descriptors["sheet_fraction_seq"] = np.nan

    # Amino acid composition
    aa_percent = analysis.amino_acids_percent

    for aa, value in aa_percent.items():
        descriptors[f"aa_comp_{aa}"] = value

    # Dipeptide composition
    amino_acids = "ACDEFGHIKLMNPQRSTVWY"

    dipeptides = {}
    for aa1 in amino_acids:
        for aa2 in amino_acids:
            dipeptides[aa1 + aa2] = 0

    for i in range(len(sequence) - 1):
        dp = sequence[i:i+2]
        if dp in dipeptides:
            dipeptides[dp] += 1

    total_dp = max(1, len(sequence) - 1)

    for dp, count in dipeptides.items():
        descriptors[f"dipep_{dp}"] = count / total_dp

    descriptors.update(residue_statistics(sequence))

    return descriptors


# =====================================================================================
# STRUCTURAL DESCRIPTORS
# =====================================================================================

def calculate_structure_descriptors(pdb_file):

    descriptors = {}

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_file)

    atoms = list(structure.get_atoms())
    residues = list(structure.get_residues())
    chains = list(structure.get_chains())

    coords = np.array([atom.coord for atom in atoms])

    descriptors["num_atoms"] = len(atoms)
    descriptors["num_residues"] = len(residues)
    descriptors["num_chains"] = len(chains)

    # Radius of gyration
    descriptors["radius_of_gyration"] = compute_radius_of_gyration(coords)

    # Bounding box
    bbox = compute_bounding_box(coords)

    descriptors["bbox_x"] = bbox[0]
    descriptors["bbox_y"] = bbox[1]
    descriptors["bbox_z"] = bbox[2]

    # Centroid
    centroid = np.mean(coords, axis=0)

    descriptors["centroid_x"] = centroid[0]
    descriptors["centroid_y"] = centroid[1]
    descriptors["centroid_z"] = centroid[2]

    # Approx volume
    descriptors["approx_volume"] = bbox[0] * bbox[1] * bbox[2]

    # FreeSASA
    try:
        fs_structure = freesasa.Structure(pdb_file)
        result = freesasa.calc(fs_structure)

        total_area = result.totalArea()

        descriptors["sasa_total"] = total_area

        try:
            total_volume = result.totalVolume()
        except:
            total_volume = np.nan

        descriptors["protein_volume"] = total_volume

        descriptors["compactness"] = compute_compactness(
            total_volume if total_volume else 1,
            total_area
        )

    except Exception as e:
        logging.warning(f"FreeSASA failed for {pdb_file}: {e}")

    # Secondary structure via DSSP
    try:
        model = structure[0]

        dssp = DSSP(model, pdb_file, dssp='mkdssp')

        ss = [dssp[key][2] for key in dssp.keys()]

        helix = sum(s in ["H", "G", "I"] for s in ss)
        sheet = sum(s in ["E", "B"] for s in ss)
        coil = sum(s == "-" for s in ss)

        total = len(ss)

        descriptors["alpha_helix_fraction"] = safe_div(helix, total)
        descriptors["beta_sheet_fraction"] = safe_div(sheet, total)
        descriptors["coil_fraction"] = safe_div(coil, total)

    except Exception as e:
        logging.warning(f"DSSP failed for {pdb_file}: {e}")

    return descriptors, structure


# =====================================================================================
# POCKET DESCRIPTORS USING FPOCKET
# =====================================================================================

def run_fpocket(pdb_file):

    pdb_path = Path(pdb_file)
    output_dir = pdb_path.with_suffix("")
    output_dir = str(output_dir) + "_out"

    # Remove old directory
    if os.path.exists(output_dir):
        shutil.rmtree(output_dir)

    cmd = [FPOCKET_EXECUTABLE, "-f", pdb_file]

    try:
        subprocess.run(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=True
        )

        return output_dir

    except Exception as e:
        logging.warning(f"fpocket failed for {pdb_file}: {e}")
        return None


def parse_fpocket_descriptors(fpocket_dir):

    descriptors = {}

    info_file = os.path.join(fpocket_dir, "pockets", "pocket1_atm.pdb")

    if not os.path.exists(info_file):
        return descriptors

    # fpocket info txt
    txt_file = os.path.join(fpocket_dir, "pockets", "pocket1_info.txt")

    if not os.path.exists(txt_file):
        return descriptors

    with open(txt_file, "r") as f:
        text = f.read()

    patterns = {
        "pocket_score": r"Pocket Score\s*:\s*([-\d\.]+)",
        "druggability_score": r"Druggability Score\s*:\s*([-\d\.]+)",
        "number_alpha_spheres": r"Number of Alpha Spheres\s*:\s*([-\d\.]+)",
        "total_sasa": r"Total SASA\s*:\s*([-\d\.]+)",
        "polar_sasa": r"Polar SASA\s*:\s*([-\d\.]+)",
        "apolar_sasa": r"Apolar SASA\s*:\s*([-\d\.]+)",
        "volume": r"Volume\s*:\s*([-\d\.]+)",
        "hydrophobic_density": r"Hydrophobic Density\s*:\s*([-\d\.]+)",
        "mean_local_hydrophobic_density": r"Mean local hydrophobic density\s*:\s*([-\d\.]+)",
        "polarity_score": r"Polarity score\s*:\s*([-\d\.]+)",
        "alpha_sphere_density": r"Alpha sphere density\s*:\s*([-\d\.]+)",
        "center_x": r"Centroid alpha sphere max dist\s*:\s*([-\d\.]+)"
    }

    for key, pattern in patterns.items():
        match = re.search(pattern, text)

        if match:
            descriptors[f"pocket_{key}"] = float(match.group(1))

    return descriptors


# =====================================================================================
# MAIN PROCESSING FUNCTION
# =====================================================================================

def process_protein(pdb_file):

    protein_name = Path(pdb_file).stem

    try:

        logging.info(f"Processing {protein_name}")

        # Structure descriptors
        structure_desc, structure = calculate_structure_descriptors(pdb_file)

        # Sequence extraction
        sequence = get_sequence_from_structure(structure)

        # Sequence descriptors
        seq_desc = calculate_sequence_descriptors(sequence)

        # Pocket descriptors
        fpocket_output = run_fpocket(pdb_file)

        if fpocket_output:
            pocket_desc = parse_fpocket_descriptors(fpocket_output)
        else:
            pocket_desc = {}

        final_desc = {
            "protein": protein_name
        }

        final_desc.update(seq_desc)
        final_desc.update(structure_desc)
        final_desc.update(pocket_desc)

        return final_desc, None

    except Exception as e:

        logging.error(f"FAILED: {protein_name}")
        logging.error(traceback.format_exc())

        return None, {
            "protein": protein_name,
            "error": str(e)
        }


# =====================================================================================
# MAIN
# =====================================================================================

def main():

    start_time = time.time()

    pdb_files = glob.glob(os.path.join(INPUT_DIR, "*.pdb"))

    logging.info(f"Found {len(pdb_files)} proteins")

    if len(pdb_files) == 0:
        print("No PDB files found.")
        return

    results = []
    failed = []

    # Multiprocessing
    if USE_MULTIPROCESSING:

        with Pool(N_CORES) as pool:

            for result, fail in tqdm(
                pool.imap_unordered(process_protein, pdb_files),
                total=len(pdb_files),
                desc="Processing Proteins"
            ):

                if result:
                    results.append(result)

                if fail:
                    failed.append(fail)

    else:

        for pdb_file in tqdm(pdb_files):

            result, fail = process_protein(pdb_file)

            if result:
                results.append(result)

            if fail:
                failed.append(fail)

    # Save results
    df = pd.DataFrame(results)

    df.to_csv(OUTPUT_CSV, index=False)

    logging.info(f"Saved descriptors to {OUTPUT_CSV}")

    # Save failures
    if len(failed) > 0:

        failed_df = pd.DataFrame(failed)

        failed_df.to_csv(FAILED_CSV, index=False)

        logging.info(f"Saved failed proteins to {FAILED_CSV}")

    runtime = time.time() - start_time

    logging.info(f"Total runtime: {runtime:.2f} seconds")
    logging.info("Processing completed successfully.")

    print("\n======================================")
    print("PROCESSING COMPLETED")
    print("======================================")
    print(f"Proteins processed : {len(results)}")
    print(f"Failed proteins    : {len(failed)}")
    print(f"Output CSV         : {OUTPUT_CSV}")
    print(f"Log file           : {LOG_FILE}")
    print("======================================\n")


# =====================================================================================
# ENTRY POINT
# =====================================================================================

if __name__ == "__main__":
    main()
