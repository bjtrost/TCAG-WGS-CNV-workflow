# TCAG-WGS-CNV-workflow
Scripts involved in our workflow for detecting CNVs from WGS data using read depth-based methods
Paper currently under review in American Journal of Human Genetics

This README file lists, and explains the purpose of, each script.
The scripts are divided into three categories:

1) "main scripts", which are designed to be called directly;
2) "accessory scripts", which are only meant to be called by/used by the main scripts; and
3) "commands", which are not meant to be called as scripts, but rather contain commands
(for running, e.g., BWA, GATK, the CNV detection algorithms) that the user
can copy-and-paste, replace placeholder filenames with their own,
and then execute sequentially on their system/cluster.

For instructions on each script, as well as example usage, please refer to the comments at the beginning
of each one.

## Main scripts (designed to be called directly)

* benchmark_overlap_counts.py: Output counts summarizing how CNVs in different categories overlap with the benchmark.
* CNV_overlap.py: Finding overlapping CNV calls from either the CNV-detection algorithms, or different benchmark methods, or both.
* CNV_read_depth_checker.sh: Calculate the ratio between the read depth of a CNV and the read depth of the same-size surrounding regions.
* compare_CNVs_to_benchmark.py: Compare CNVs output from the CNV-detection algorithms to a CNV benchmark.
* compare_with_RLCR_definition.py: Compare CNV calls that have been converted to the common format with the RLCR definition. Requires the "intervaltree" Python module to be installed.
* convert_CNV_calls_to_common_format.py: Convert CNV calls to common format.
* IQR_samtools_depth.sh: Calculates IQR from a BAM file.
* merge_Genome_STRiP.py: Use this script on a Genome STRiP file that has already been converted to the common format using convert_CNV_calls_to_common_format.py in order to merge overlapping calls.
* process_cnvs.erds+.sh: Use this script to perform the CNVnator-ERDS merging that was used in Stage 3 of the study (the analysis of rare, genic CNVs in the Autism Speaks MSSNG cohort).
* reproduce_results.sh: A script that runs all the other main scripts in an appropriate sequence.
* split_HuRef_benchmark.sh: Split the file containing the HuRef benchmark CNVs into separate files, one for each benchmark technology.

## Accessory scripts (NOT to be called directly)
* add_features.py: Used by process_cnvs.erds+.sh.
* Canvas.py: Custom library for converting Canvas output to common format.
* cnMOPS.py: Custom library for converting cn.MOPS output to common format.
* CNVnator.py: Custom library for converting CNVnator output to common format.
* CNVworkflowlib.py: Custom library of python functions used by other python scripts.
* ERDS.py: Custom library for converting ERDS output to common format.
* format_cnvnator_results.py: Used by process_cnvs.erds+.sh.
* format_erds_results.py: Used by process_cnvs.erds+.sh.
* functions.py: Used by process_cnvs.erds+.sh.
* Genome_STRiP.py: Custom library for converting Genome_STRiP output to common format.
* get_normalized_depth.py: Used by CNV_read_depth_checker.sh to actually calculate the normalized depth of a CNV.
* index_samtools_depth.py: Used by CNV_read_depth_checker.sh to index a "samtools depth" file for fast use of get_normalized_depth.py.
* IQR_samtools_depth.py: Does most of the work involved in calculating IQR from a BAM file.
* merge_cnvnator_results.py: Used by process_cnvs.erds+.sh.
* merge_erds_results.py: Used by process_cnvs.erds+.sh.
* myvcf.py: Custom library of python functions for dealing with VCF files.
* RDXplorer.py: Custom library for converting RDXplorer output to common format.
* SV.py: Custom library of python functions for representing SVs/CNVs.

## Commands (designed for the user to execute commands one-by-one)
* commands.sh: a list of commands for running BWA, GATK, etc.
