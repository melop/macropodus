TRINITY=/data/software/trinityrnaseq-v2.12.0/Trinity

$TRINITY --CPU 64 --seqType fq --max_memory 100G --samples_file samples.txt --jaccard_clip --output trinity_assm  > out.log 2>&1 
