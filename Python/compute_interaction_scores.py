import os
import argparse
import logging
import numpy as np
import pandas as pd
import retrieveLigRecInfo as LRinfo
import RaCInG_input_generation as input_gen


def compute_interaction_scores(weight_type, cellfrac_path, lr_weights_dir, lr_pairs_dir, source_cell_type="Tumor", target_cell_type="M", output_dir="", cancer_name="GBM", is_directed = True):
    liglist, reclist, Dcell, Dconn, cellnames, ligands, receptors, _ = input_gen.generateInput(
        weight_type, cancer_name, cellfrac_path=cellfrac_path, lr_weights_dir=lr_weights_dir, lr_pairs_dir=lr_pairs_dir)

    # Compute interaction scores
    if (is_directed):
        weights = calculateBaysianInteractionScoreDirected(liglist, reclist, Dcell, Dconn, cellnames, source_cell_type, target_cell_type)

    else:
        weights = LRinfo.calculateBaysianInteractionScoreUndirected(
        liglist, reclist, Dcell, Dconn, cellnames, source_cell_type, target_cell_type)



    # Save all results
    np.save(
        file=f"{output_dir}/{source_cell_type}__{target_cell_type}__weights.npz", arr=weights)

    np.save(
        file=f"{output_dir}/{source_cell_type}__{target_cell_type}__ligands.npz", arr=ligands)

    np.save(
        file=f"{output_dir}/{source_cell_type}__{target_cell_type}__receptors.npz", arr=receptors)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compute Bayesian scores")
    parser.add_argument("-o", "--output_dir", type=str,
                        help="Path to output files")
    parser.add_argument("-c", "--cancer_type", type=str,
                        help="Cancer type", default="GBM")
    parser.add_argument("-lrp", "--lr_pairs_dir", type=str,
                        help="Path to LR pairs directory", default="")
    parser.add_argument("-cf", "--cellfrac_path", type=str,
                        help="Path to cell fraction file", default="")
    parser.add_argument("-lrw", "--lr_weights_dir", type=str,
                        help="Path to LR dir weights directory", default="")
    parser.add_argument("-w", "--weight_type", type=str,
                        help="Weight type", default="min")
    parser.add_argument("-s", "--source_cell_type", type=str,
                        default="Tumor", help="Source cell type")
    parser.add_argument("-t", "--target_cell_type", type=str,
                        default="M", help="Target cell type")
    args = parser.parse_args()

    # Setup logging
    logging.basicConfig(level=logging.DEBUG, format='%(asctime)s %(message)s')

    if (not os.path.exists(args.output_dir)):
        os.makedirs(args.output_dir)

    compute_interaction_scores(
        args.weight_type, args.cancer_type, args.cellfrac_path, args.lr_weights_dir,)
