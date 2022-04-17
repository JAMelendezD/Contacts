
mkdir $4
frames=$(python get_frames.py $1 $2 2>&1)  

#sele1="bynum 1:6557"
#sele2="bynum 13115:14918"
sele1="bynum 6558:13114"
sele2="bynum 14919:16720"
startres1=184
resbefore2=816
dt=0.02
cutoff=0.1
chain1=A
chain2=C

python contacts.py $1 $2 0 999999 5.0 "$sele1" "$sele2" ./${4}/phobic >> ${4}.log 2>/dev/null
python salt_bridges.py $1 $2 0 999999 4.0 "$sele1" "$sele2" ./${4}/salt >> ${4}.log 2>/dev/null
python hbonds.py $1 $2 0 999999 3.5 30.0 "$sele1" "$sele2" ./${4}/hbonds >> ${4}.log 2>/dev/null

for name in phobic salt hbonds
do
    for i in 0 1
    do
       python count.py ./${4}/${name}.txt $startres1 -$resbefore2 $frames ./${4}/${name}_count --mode $i
    done
    if [ $name == hbonds ]
    then
        python time_series.py ./${4}/${name}.txt ./${4}/${name}_count.dat $startres1 -$resbefore2 $frames $dt $cutoff ./${4}/${name}_series
    else
        python time_series.py ./${4}/${name}.txt ./${4}/${name}_count.dat $startres1 -$resbefore2 $frames $dt $cutoff --mode 1 ./${4}/${name}_series
    fi
done

python network.py ./${4}/hbonds_count.dat ./${4}/hbond_salt 0.1 --auxmark 1 --auxfile ./${4}/salt_count.dat --rowsep 0.7 >> ${4}.log 2>/dev/null
python network.py ./${4}/phobic_count.dat ./${4}/phobic 0.1 --rowsep 0.7 >> ${4}.log 2>/dev/null

python paint_pdb.py $3 ./${4}/hbonds_count_1.dat $chain1 ./${4}/tmp.pdb >> ${4}.log 2>/dev/null
python paint_pdb.py ./${4}/tmp.pdb ./${4}/hbonds_count_2.dat $chain2 ./${4}/final.pdb >> ${4}.log 2>/dev/null


rm ./${4}/*.dot
rm ./${4}/tmp.pdb
