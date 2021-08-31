import filecmp
import os
import argparse

def compare_with_gold_file(filename, out_dir, gold_dir):
    out_filename  = os.path.join(out_dir, filename)
    gold_filename = os.path.join(gold_dir, "gold_{:s}".format(filename))
    if filecmp.cmp(out_filename, gold_filename, shallow = False):
        print (f"OK {filename}")
    else:
        print (f"Error in {filename}")


def parse_args():
    parser = argparse.ArgumentParser(description='This code compares Tejaas example output with the known result.')
    parser.add_argument('--outdir',
                        type=str,
                        dest='outdir',
                        metavar='DIR',
                        default='.',
                        help='Output directory of Tejaas example')
    return parser.parse_args()



if __name__ == "__main__":
    opts     = parse_args()
    cur_dir  = os.getcwd()
    gold_dir = os.path.join(cur_dir, "gold")
    out_dir  = os.path.join(os.path.realpath(opts.outdir), "data")
    for filename in ["result_rr.txt", 
                     "result_gene_snp_list.txt", 
                     "result_gene_snp_list_knn.txt"]:
        compare_with_gold_file(filename, out_dir, gold_dir)
