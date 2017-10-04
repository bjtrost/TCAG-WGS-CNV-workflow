#!/usr/bin/env bash

# Convert raw output files to common format
./convert_CNV_calls_to_common_format.py test_data/89931_500.HAS.Canvas.vcf Canvas > Canvas.txt
./convert_CNV_calls_to_common_format.py test_data/89931_500.cn.MOPS.txt cn.MOPS > cn.MOPS.txt
./convert_CNV_calls_to_common_format.py test_data/89931_500.CNVnator.bin500.txt CNVnator > CNVnator.txt
./convert_CNV_calls_to_common_format.py test_data/89931_500.erds.del200.vcf ERDS > ERDS.txt
./convert_CNV_calls_to_common_format.py test_data/89931_500.GenomeSTRiP.CNVpipeline.vcf Genome_STRiP > Genome_STRiP.tmp
./convert_CNV_calls_to_common_format.py test_data/89931_500.RDX.txt RDXplorer > RDXplorer.txt

# Merge overlapping Genome STRiP calls
./merge_Genome_STRiP.py Genome_STRiP.tmp > Genome_STRiP.txt
rm Genome_STRiP.tmp

# Compare each algorithm against HuRef benchmark
for algorithm in Canvas cn.MOPS CNVnator ERDS Genome_STRiP RDXplorer; do
    ./compare_CNVs_to_benchmark.py $algorithm.txt ../supplemental_files/File_S1.txt > $algorithm.txt.bm
done

# Compare benchmark-compared files against RLCR definition
./compare_with_RLCR_definition.py test_data/RLCRs.txt Canvas.txt.bm cn.MOPS.txt.bm CNVnator.txt.bm ERDS.txt.bm Genome_STRiP.txt.bm RDXplorer.txt.bm

# Make .counts files, which give different categories of overlap with repeat definition and RLCR definition
for file in *.RLCR; do
    ./benchmark_overlap_counts.py $file > $file.counts

# Make file containing overlap of all methods
./CNV_overlap.py Canvas.txt,cn.MOPS.txt,CNVnator.txt,ERDS.txt,Genome_STRiP.txt,RDXplorer.txt File_S1.txt.Affymetrix6.0,File_S1.txt.Agilent24M,File_S1.txt.Illumina1M,File_S1.txt.NimbleGen42M,File_S1.txt.SplitRead,File_S1.txt.assembly.comparison,File_S1.txt.cg.cnv,File_S1.txt.cg.mp > all_overlap.txt
