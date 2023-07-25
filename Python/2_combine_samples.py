import execute_sim as exSim
import argparse
import logging

parser = argparse.ArgumentParser(description="Combine simulations from different samples")

parser.add_argument("-f", "--input_dir", type=str, help="Path to input files")
parser.add_argument("-o", "--output_dir", type=str,
                    help="Path to output files")
parser.add_argument("-c", "--cancer_type", type=str, help="Cancer type")
parser.add_argument("-av", "--av", type=int,
                    help="Average degree of network to be generated", default=2
                    )
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

    exSim.combinePatients(input_dir = args.input_dir, cancer_type = args.cancer_type, communication_type = args.communication_type, norm = args.norm_dist, av = args.av, output_dir = args.output_dir)