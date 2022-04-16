# python_tests


```
python contacts.py protein.tpr protein_nj_fitbb.xtc 0 999999 5.0 ./results/phobic_interface1_run1.txt
python salt.py protein.tpr protein_nj_fitbb.xtc 0 999999 4.0 ./results/salt_interface1_run1.txt
python hbonds.py protein.tpr protein_nj_fitbb.xtc 0 999999 3.5 30.0 ./results/hbonds_interface1_run1.txt
```

```
python count.py ./results/hbonds_interface1_run1.txt 184 -816 26102 ./results/hbonds_count_interface1_run1 --mode 0
python count.py ./results/hbonds_interface1_run1.txt 184 -816 26102 ./results/hbonds_count_interface1_run1 --mode 1
```

```
python time_series.py ./results/hbonds_interface1_run1.txt 26102 20 5.0 ./results/series_hbonds
```

```
python network.py ./results/hbonds_count_interface1_run1.dat hbond_salt 0.1 --auxmark 1 --auxfile ./results/salt_count_interface1_run1.dat --rowsep 0.7
python network.py ./results/phobic_count_interface1_run1.dat phobic 0.1 --rowsep 0.7
```

```
python paint_pdb.py protein.pdb hbonds_count_interface1_run1_1.dat A tmp.pdb
python paint_pdb.py tmp.pdb hbonds_count_interface1_run1_2.dat B final.pdb
```