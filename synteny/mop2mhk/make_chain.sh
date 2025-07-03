#mop coordinate to mhk coordinate, mop is reference genome and mhk is query genome in chain file

query=/data/projects/rcui/mhk/nextdenovo/nextpolish_stlfr2/mhk.LG.nextdenovo.sgs.2.sorted.fa
ref=/data/projects/rcui/mop/nextpolish/nextpolish_stlfr2/nextpolish/mop.LG.nextdenovo.sgs.2.fa
transanno=/data/software/transanno/transanno-x86_64-unknown-linux-musl-v0.4.4/transanno

minimap2 -cx asm5 --cs $query $ref > mop2mhk.paf

$transanno minimap2chain mop2mhk.paf --output mop2mhk.chain
