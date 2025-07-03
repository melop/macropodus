#caculate ave dxy of whole genome
cat pixy_output_dxy.txt | cut -f1,2 |sort|uniq > pops_dxy.list
cat pops_dxy.list | while read line
do
    a=$(echo $line | awk '{print $1}');
    b=$(echo $line | awk '{print $2}');
    cat pixy_output_dxy.txt |awk '$1 == "'$a'" && $2 == "'$b'" {suma+=$8;sumb+=$9} END {print "'$a'", "'$b'", suma, sumb, suma/sumb}' >> dxy_average.txt
done

#caculate pi of whole genome
cat pixy_output_pi.txt |cut -f1 |sort|uniq > pops_pi.list
cat pops_pi.list | while read line
do
    a=$(echo $line | awk '{print $1}');
    cat pixy_output_pi.txt | awk '$1 == "'$a'" {suma+=$7;sumb+=$8} END {print "'$a'", suma, sumb, suma/sumb}' >> pi_average.txt
done
