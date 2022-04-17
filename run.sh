
startres1=184
resbefore2=816
frames=26102
dt=0.02
cutoff=0.1 

python contacts.py $1 $2 0 999999 5.0 ./${4}/phobic_${5}
python salt_bridges.py $1 $2 0 999999 4.0 ./${4}/salt_${5}
python hbonds.py $1 $2 0 999999 3.5 30.0 ./${4}/hbonds_${5} 

for name in phobic salt hbonds
do
    for i in 0 1
    do
       python count.py ./${4}/${name}_${3}.txt $startres1 -$resbefore2 $frames ./${4}/${name}_count_${5} --mode $i
    done
    if [ $name == hbonds ]
    then
        python time_series.py ./${4}/${name}_${5}.txt ./${4}/${name}_count_${5}.dat $startres1 -$resbefore2 $frames $dt $cutoff ./${4}/${name}_series
    else
        python time_series.py ./${4}/${name}_${5}.txt ./${4}/${name}_count_${5}.dat $startres1 -$resbefore2 $frames $dt $cutoff --mode 1 ./${4}/${name}_series
    fi
done

python network.py ./${4}/hbonds_count_${5}.dat ./${4}/hbond_salt 0.1 --auxmark 1 --auxfile ./${4}/salt_count_${5}.dat --rowsep 0.7
python network.py ./${4}/phobic_count_${5}.dat ./${4}/phobic 0.1 --rowsep 0.7

python paint_pdb.py $3 ./${4}/hbonds_count_${5}_1.dat A ./${4}/tmp.pdb
python paint_pdb.py ./${4}/tmp.pdb ./${4}/hbonds_count_${5}_2.dat C ./${4}/final_${5}.pdb
