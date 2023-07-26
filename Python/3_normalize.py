import os
import argparse
import logging

import txt_to_csv as convert

parser = argparse.ArgumentParser(description="Generate Network")

parser.add_argument("-f", "--input_dir", type=str, help="Path to input files")
parser.add_argument("-o", "--output_dir", type=str,
                    help="Path to output files")
parser.add_argument("-c", "--cancer_type", type=str, help="Cancer type")
parser.add_argument("-p", "--patient_id", type=int,
                    help="Id of patient (number)", default=1)

parser.add_argument("-n", "--n_cells", type=int,
                    help="Number of cells in network to be generated", default=100)

parser.add_argument("-av", "--av", type=int,
                    help="Average degree of network to be generated", default=2
                    )
parser.add_argument("-itNo", "--n_graphs", type=int,
                    help="Number of graphs to be generated", default=10)
parser.add_argument("-w", "--weight_type", type=str,
                    help="Weight type", default="min")

parser.add_argument("-ct", "--communication_type", type=str,
                    help="Communication type: direct(D), wedge(W), trust triangle(TT), or cycle triangle(CT)", default="W")
parser.add_argument("-norm", "--norm_dist", type=int,
                    help="Normalise interaction distribution", default=0)
args = parser.parse_args()


if __name__ == "__main__":
    # Setup logging
    logging.basicConfig(level=logging.DEBUG, format='%(asctime)s %(message)s')

    if (not os.path.exists(args.output_dir)):
        os.makedirs(args.output_dir)

    convert.Generate_normalised_count_csv(
        cancer_type=args.cancer_type, weight_type=args.weight_type, triangle_types=args.communication_type, average=args.av, noCells=args.n_cells, input_dir=args.input_dir, remove_direction=True, output_dir=args.output_dir)
