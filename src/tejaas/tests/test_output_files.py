import unittest
import os
import subprocess
import pathlib

import tejaas.tester
  
class TestOutputFiles(unittest.TestCase):


    def __init__(self, *args, **kwargs):
        super(TestOutputFiles, self).__init__(*args, **kwargs)
        self.test_dir = pathlib.Path(__file__).parent.resolve()
        self.out_dir  = os.path.join(self.test_dir, "data")
        self.gold_dir = os.path.join(self.test_dir, "gold")
        self._run_bash_script(os.path.join(self.test_dir, "run_mwe.sh"))


    def test_rr_outfiles(self):
        for filename in ["result_rr.txt", "result_gene_snp_list.txt", "result_gene_snp_list_knn.txt"]:
            self._compare_with_gold_file(filename)
        return


    #def test_jpanull(self):
    #    self._compare_with_gold_file("result_jpanull.txt")
    #    return


    #def test_jpa_outfiles(self):
    #    for filename in ["result_fr_jpa.txt", "result_fr_gene_snp_list.txt"]:
    #        self._compare_with_gold_file(filename)
    #    return


    def _compare_with_gold_file(self, filename):
        out_filename = os.path.join(self.out_dir, filename)
        gold_filename = os.path.join(self.gold_dir, "gold_{:s}".format(filename))
        self._compare_files(out_filename, gold_filename)
        return
            

    def _compare_files(self, file1, file2):
        f1 = open(file1, "r")
        f2 = open(file2, "r")
        self.assertMultiLineEqual(f1.read(), f2.read())
        f1.close()
        f2.close()
        return


    def _run_bash_script(self, script_file):
        cmd     = ["bash", script_file, self.test_dir]
        with subprocess.Popen(cmd, stdout = subprocess.PIPE, stderr = subprocess.PIPE, bufsize=1, universal_newlines=True) as process:
            for line in process.stdout:
                print (line, end = "")
        retcode = process.returncode
        self.assertEqual(retcode, 0)
        return        

if __name__ == '__main__':
    tejaas.tester.main()
