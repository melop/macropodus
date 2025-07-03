zcat /data/projects/shareddata/mhk.stLFR/raw_data_2021-02-02/*/barcode*/qc_report/fastq/*.pass*gz | pigz -p 64 -c > mhk.pooled.ont.isoseq.fastq.gz
