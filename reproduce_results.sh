#!/usr/bin/env bash

# Directory $test_data_dir should have contents of supplementary files (File_S1.txt, *.HuRef.*, and RLCRs.txt)

output_dir="reproduce_results_out"
mkdir $output_dir

test_data_dir="test_data"

# Convert raw output files to common format
./convert_CNV_calls_to_common_format.py $test_data_dir/Canvas.HuRef.vcf Canvas > $output_dir/Canvas.txt
./convert_CNV_calls_to_common_format.py $test_data_dir/cn.MOPS.HuRef.txt cn.MOPS > $output_dir/cn.MOPS.txt
./convert_CNV_calls_to_common_format.py $test_data_dir/CNVnator.HuRef.txt CNVnator > $output_dir/CNVnator.txt
./convert_CNV_calls_to_common_format.py $test_data_dir/ERDS.HuRef.vcf ERDS > $output_dir/ERDS.txt
./convert_CNV_calls_to_common_format.py $test_data_dir/Genome_STRiP.HuRef.vcf Genome_STRiP > $output_dir/Genome_STRiP.tmp
./convert_CNV_calls_to_common_format.py $test_data_dir/RDXplorer.HuRef.txt RDXplorer > $output_dir/RDXplorer.txt

# Merge overlapping Genome STRiP calls
./merge_Genome_STRiP.py $output_dir/Genome_STRiP.tmp > $output_dir/Genome_STRiP.txt
rm $output_dir/Genome_STRiP.tmp

# Compare each algorithm against HuRef benchmark
for algorithm in Canvas cn.MOPS CNVnator ERDS Genome_STRiP RDXplorer; do
    ./compare_CNVs_to_benchmark.py $output_dir/$algorithm.txt $test_data_dir/File_S1.txt > $output_dir/$algorithm.txt.bm
done

# Compare benchmark-compared files against RLCR definition
./compare_with_RLCR_definition.py $test_data_dir/RLCRs.txt $output_dir/Canvas.txt.bm $output_dir/cn.MOPS.txt.bm $output_dir/CNVnator.txt.bm $output_dir/ERDS.txt.bm $output_dir/Genome_STRiP.txt.bm $output_dir/RDXplorer.txt.bm

# Make .counts files, which give different categories of overlap with repeat definition and RLCR definition
for file in *.RLCR; do
    ./benchmark_overlap_counts.py $output_dir/$file > $output_dir/$file.counts

# Make files containing CNVs from individual benchmark methods
cp $test_data_dir/File_S1.txt $output_dir
./split_HuRef_benchmark.sh $output_dir/File_S1.txt

# Make file containing overlap of all methods
cd $output_dir
../CNV_overlap.py Canvas.txt,cn.MOPS.txt,CNVnator.txt,ERDS.txt,Genome_STRiP.txt,RDXplorer.txt File_S1.txt.Affymetrix6.0,File_S1.txt.Agilent24M,File_S1.txt.Illumina1M,File_S1.txt.NimbleGen42M,File_S1.txt.SplitRead,File_S1.txt.assembly.comparison,File_S1.txt.cg.cnv,File_S1.txt.cg.mp > all_overlap.txt
cd ..
