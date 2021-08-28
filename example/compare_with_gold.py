import filecmp
import os

def compare_with_gold_file(filename, out_dir, gold_dir):
    out_filename  = os.path.join(out_dir, filename)
    gold_filename = os.path.join(gold_dir, "gold_{:s}".format(filename))
    if filecmp.cmp(out_filename, gold_filename, shallow = False):
        print (f"OK {filename}")
    else:
        print (f"Error in {filename}")


if __name__ == "__main__":
    cur_dir = os.getcwd()
    out_dir = os.path.join(cur_dir, "data")
    gold_dir = os.path.join(cur_dir, "gold")
    for filename in ["result_rr.txt", 
                     "result_gene_snp_list.txt", 
                     "result_gene_snp_list_knn.txt"]:
        compare_with_gold_file(filename, out_dir, gold_dir)
