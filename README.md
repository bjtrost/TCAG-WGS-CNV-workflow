# TCAG-WGS-CNV-workflow
Scripts involved in our workflow for detecting CNVs from WGS data using read depth-based methods

This file lists, and explains the purpose of, each script in this directory.
The scripts are divided into three categories:

1) "main scripts", which are designed to be called directly;
2) "accessory scripts", which are only meant to be called by/used by the main scripts; and
3) "commands", which are not meant to be called as scripts, but rather contain commands
(e.g., BWA, GATK) that the user can copy-and-paste, replace placeholder filenames with their own,
and then execute sequentially on their system/cluster.

For instructions on each script, as well as example usage, please refer to the comments at the beginning
of each one. If you have any difficulties setting up and running these scripts, or if you encounter any bugs,
please e-mail brett.trost@sickkids.ca.

# Main scripts (designed to be called directly)
* IQR_samtools_depth.sh: Calculates IQR from a BAM file
* convert_CNV_calls_to_common_format.py: Convert CNV calls to common format
* merge_Genome_STRiP.py: Use this script on a Genome STRiP file that has already been converted to the common format using convert_CNV_calls_to_common_format.py in order to merge overlapping calls.
* compare_with_RLCR_definition.py: Compare CNV calls that have been converted to the common format with the RLCR definition. Requires the "intervaltree" Python module to be installed.

# Accessory scripts (NOT to be called directly)
* CNVworkflowlib.py: Custom library of python functions used by other python scripts.
* myvcf.py: Custom library of python functions for dealing with VCF files
* SV.py: Custom library of python functions for representing SVs/CNVs
* IQR_samtools_depth.py: Does most of the work involved in calculating IQR from a BAM file

* Canvas.py: Custom library for converting Canvas output to common format
* cnMOPS.py: Custom library for converting cn.MOPS output to common format
* CNVnator.py: Custom library for converting CNVnator output to common format
* ERDS.py: Custom library for converting ERDS output to common format
* Genome_STRiP.py: Custom library for converting Genome_STRiP output to common format
* RDXplorer.py: Custom library for converting RDXplorer output to common format

# Commands (designed for the user to execute commands one-by-one)
* commands.sh: a list of commands for running BWA, GATK, etc.
