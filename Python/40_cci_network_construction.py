# Estimated running time: 1 second
import scprep
import pandas as pd
import cell2cell as c2c
import logging
import argparse
import os
import network_generation as netGen
import RaCInG_input_generation as readIn
import numpy as np


# Import modules

parser = argparse.ArgumentParser(
    description="Generate Network")

# Parse arguments
parser.add_argument("-f", "--input_dir", type=str,
                    help="Path to input files")
parser.add_argument("-o", "--output_dir", type=str,
                    help="Path to output files")
parser.add_argument("-p", "--patient_id", type=int,
                    help="Id of patient (number)")
parser.add_argument("-N", "--N", type=int,
                    help="Number of cells in network to be generated")
parser.add_argument("-avdeg", "--avdeg", type=int,
                    help="Average degree of network to be generated")
args = parser.parse_args()

# Making the LR-distribution uniform


def convert_to_uniform(LRdistr):
    LRdistrUniform = np.zeros_like(LRdistr)
    normvec = 1 / np.count_nonzero(LRdistr, axis=(0, 1))
    for i in range(LRdistr.shape[2]):
        copy = LRdistr[:, :, i].copy()
        copy[copy > 0] = normvec[i]
        LRdistrUniform[:, :, i] = copy
    return LRdistrUniform


def main():
    np.random.seed(args.seed)
    # Create output directories if they don't exist
    # if (args.output_dir == "output"):
    #     args.output_dir = f"{args.work_dir}/output"
    if (args.N > 1e5):
        output_dir = f"{args.output_dir}/{args.cancer_type}/{args.N:.2e}_{args.avdeg}_{'norm' if args.is_uniform_lr_dist else ''}"
    else:
        output_dir = f"{args.output_dir}/{args.cancer_type}/{args.N}_{args.avdeg}_{'norm' if args.is_uniform_lr_dist else ''}"
    if (not os.path.exists(output_dir)):
        os.makedirs(output_dir)

    # Setup logging
    logging.basicConfig(level=logging.DEBUG, format='%(asctime)s %(message)s')

    logging.info("Generate input...")
    Lmatrix, Rmatrix, Cdistr, LRdistr, cellTypes, ligs, recs, _ = readIn.generateInput(
        "min", "BLCA", folder_dir=args.input_dir)
    np.random.seed(1)

    if (args.is_uniform_lr_dist):
        LRdistr = convert_to_uniform(LRdistr)

    # Note we use only data of patient 0 from the input data
    V, E, types = netGen.model1(args.N, args.avdeg, Lmatrix, Rmatrix,
                                Cdistr[args.patient_id, :], LRdistr[:, :, args.patient_id], genRandom=False)
    # Save network
    np.savetxt(
        f"{output_dir}/vertices.csv", V, delimiter=",")

    np.savetxt(f"{output_dir}/edges.csv",
               E, delimiter=",")
    np.savetxt(f"{output_dir}/lr_info_edges.csv",
               types, delimiter=",")
    return V, E, types


if __name__ == "__main__":
    args.input_dir = "/Users/joankant/Library/CloudStorage/OneDrive-UHN/RaCInG/output/BLCA"
    args.output_dir = "/Users/joankant/Library/CloudStorage/OneDrive-UHN/RaCInG/output/40_cci_network_construction"

    # Patient for whom to generate the network
    args.patient_id = 0
    # Network properties:
    # N = number of cells
    # avdeg = average degree
    args.N = 20
    args.avdeg = 2
    args.seed = 1
    args.is_uniform_lr_dist = True
    args.cancer_type = "BLCA"

    V, E, types = main()
