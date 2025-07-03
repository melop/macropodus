genomeSize=438000000
#singularity run docker://dfam/tetools:latest calcDivergenceFromAlign.pl ../scf.fa.cat.gz  -s summary.txt

singularity run docker://dfam/tetools:latest  createRepeatLandscape.pl -div summary.txt -g $genomeSize | sed 's|www.google.com/jsapi|www.gstatic.com/charts/loader.js|' > mop.landscape.html
