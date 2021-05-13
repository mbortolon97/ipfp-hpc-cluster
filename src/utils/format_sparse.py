from scipy.sparse import csr_matrix, load_npz
import argparse

def input():
    parser = argparse.ArgumentParser(description="Convert a sparse SciPy matrix into the format for the IPFP")
    parser.add_argument("input_filepath", type=str, help="the filepath where the SciPy ipfp matrix are stored")
    parser.add_argument("output_filepath", type=str, help="the filepath where the IPFP matrix will be saved")
    args = parser.parse_args()
    input_filepath = args.input_filepath
    output_filepath = args.output_filepath

    return input_filepath, output_filepath

def main(input_filepath, output_filepath):
    result = load_npz(input_filepath)
    if type(result) is csr_matrix:
        result = result.tocoo()

    
    
    with open(output_filepath, 'w') as outfile:
        outfile.write("{} {} {}\n".format(result.shape[0], result.shape[1], len(result.data)))
        for idx, data in enumerate(result.data):
            outfile.write("{} {} {}\n".format(result.row[idx], result.col[idx], data))

if __name__ == "__main__":
    input_filepath, output_filepath = input()
    main(input_filepath, output_filepath)