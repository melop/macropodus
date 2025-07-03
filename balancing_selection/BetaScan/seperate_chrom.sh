sPop=MOPdsy
sF=/data2/projects/zwang/m.hk/ROH/BetaScan/$sPop/$sPop.beta.txt.gz

less "$sF" | awk 'BEGIN { n=1; last=0; lines=""; }
{
    if ($1 > last) {
        lines = lines $0 "\n";
        last = $1;
    } else {
    	filename = "chrom_" n ".beta.txt";
	sub(/\n$/, "", lines);
	print lines > filename;
        n += 1;
        if (n > 23) {
            exit;
        }
        lines = $0 "\n";  # 重置lines为当前行
        last = $1;
    }
}' 
#>> ./test.txt

#END {
#    print n, lines;
#}


