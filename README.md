# Contact Analysis

## Generate diferent contacts from trajectory

Creates text files with the raw contact data from the trajectory

```
python contacts.py ./data/protein.tpr ./data/protein_nj_fitbb.xtc 0 999999 5.0 ./results/phobic_interface1_run1
python salt_bridges.py ./data/protein.tpr ./data/protein_nj_fitbb.xtc 0 999999 4.0 ./results/salt_interface1_run1
python hbonds.py ./data/protein.tpr ./data/protein_nj_fitbb.xtc 0 999999 3.5 30.0 ./results/hbonds_interface1_run1
```
## Count the number of contacts individually and by pairs

Using the generated text files creates a dat file with the count by pairs mode 0 or individual mode 1.

```
python count.py ./results/hbonds_interface1_run1.txt 184 -816 26102 ./results/hbonds_count_interface1_run1 --mode 0
python count.py ./results/hbonds_interface1_run1.txt 184 -816 26102 ./results/hbonds_count_interface1_run1 --mode 1
```
## Generate a time series

```
python time_series.py ./results/hbonds_interface1_run1.txt ./results/hbonds_count_interface1_run1.dat 184 -816 26102 0.02 0.1 ./results/hbond_series
```
## Create network maps

```
python network.py ./results/hbonds_count_interface1_run1.dat ./results/hbond_salt 0.1 --auxmark 1 --auxfile ./results/salt_count_interface1_run1.dat --rowsep 0.7
python network.py ./results/phobic_count_interface1_run1.dat ./results/phobic 0.1 --rowsep 0.7
```
## Map value to a pdb file

```
python paint_pdb.py ./data/protein.pdb ./results/hbonds_count_interface1_run1_1.dat A tmp.pdb
python paint_pdb.py ./results/tmp.pdb ./results/hbonds_count_interface1_run1_2.dat B final.pdb
```