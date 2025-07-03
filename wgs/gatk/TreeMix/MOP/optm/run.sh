FILE=../MOP.noMissing.treemix.frq.gz
treemix=/data/software/treemix/src/treemix
nCPU=40

nCount=0;

for m in {1..10}
   do
   for i in {1..10}
      do
      nCount=$((nCount+1));
      if [ $nCount -gt $nCPU ]; then
        wait;
        nCount=1;
      fi
      k=$((i*200));
      if [ -s test.${i}.${m}.llik ]; then
        echo test.${i}.${m}.llik done
      else
          $treemix -i $FILE -k $k -global -root MOPyn -m $m -o test.${i}.${m} > log.${i}.${m}.log 2>&1 &
      fi
      done 
done

m=0
for i in {1..10}
do
      nCount=$((nCount+1));
      if [ $nCount -gt $nCPU ]; then
        wait;
        nCount=1;
      fi
      
      k=$((i*200));
      if [ -s test.${i}.${m}.llik ]; then
        echo test.${i}.${m}.llik done
      else
          $treemix -i $FILE -k $k -global -root MOPyn -o test.${i}.${m} > log.${i}.${m}.log 2>&1 &
      fi
done

wait
