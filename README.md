    EXPRFILE="GTEx_wholeBlood_Normalzed_NoPEER_lmcorrected.txt"
    mpirun -n 4 python main.py --gx ${EXPRFILE} --simulate --test --method rr --null maf
