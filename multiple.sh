
for i in 1 2 3 4 5 6
do
   ./run.sh ./data/run${i}/protein.tpr ./data/run${i}/protein_nj_fitbb.xtc ./data/protein.pdb run${i}_interface1 &
done

