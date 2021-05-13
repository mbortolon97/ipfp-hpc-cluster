import argparse
import numpy as np
import sys

def input():
    parser = argparse.ArgumentParser(description="Convert a dense NumPy matrix into the format for the IPFP")
    parser.add_argument("input_filepath", type=str, help="the filepath where the NumPy matrix are stored")
    parser.add_argument("output_filepath", type=str, help="the filepath where the IPFP dense matrix will be saved")
    args = parser.parse_args()
    input_filepath = args.input_filepath
    output_filepath = args.output_filepath

    return input_filepath, output_filepath

def main(input_filepath, output_filepath):
    result = np.load(input_filepath)
    
    
    with open(output_filepath, 'w') as outfile:
        shape_str = ' '.join([str(x) for x in result.shape])
        outfile.write("{}\n".format(shape_str))
        for i in range(result.shape[0]):
            for j in range(result.shape[1]):
                outfile.write("{} {} {}\n".format(i, j, result[i, j]))

if __name__ == "__main__":
    input_filepath, output_filepath = input()
    main(input_filepath, output_filepath)