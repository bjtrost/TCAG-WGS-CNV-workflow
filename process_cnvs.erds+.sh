#!/usr/bin/env bash

# To start, set "python_path" to be the path to Python 2 on your system, and
# set "root" to be the directory containing this script and the other provided scripts.

# Then, make a new directory, which we will call "mycnvs". Then make the directories mycnvs/erds, mycnvs/cnvn, mycnvs/erds/original, and mycnvs/cnvn/original.
# Inside the mycnvs/erds/original directory, add files called <myid>.erds.vcf. ***Note that this is the VCF file outputted by ERDS.***
# Inside the mycnvs/cnvn/original directory, add files called <myid>.calls.txt. *** Note that this is the tab-delimited file outputted by CNVnator, not the optional VCF file ***
# For the same genome, <myid> should be the same. These tell the script which files to merge. Thus, if you have ERDS and CNVnator calls from
# two different genomes called genome1 and genome2, then mycnvs/erds/original should contain files called genome1.erds.vcf and genome2.erds.vcf,
# and mycnvs/cnvn/original should contain files called genome1.calls.txt and genome2.calls.txt. Then merging will be done for genome1.erds.vcf and genome1.calls.txt,
# as well as for genome2.erds.vcf and genome2.calls.txt.


cnv_dir=$1
python_path=/Users/btrost/scripts/python2 # Replace with actual path to Python 2
########
root=/Users/btrost/Desktop/test/TCAG-WGS-CNV-workflow/
#######
cnvn_format="$root"format_cnvnator_results.py
erds_format="$root"format_erds_results.py
erds_merge="$root"merge_erds_results.py
cnvn_merge="$root"merge_cnvnator_results.py
erds_plus="$root"add_features.py
gap_file="$root"hg19_gap.bed
echo $cnvn_format $erds_format $cnvn_merge $erds_merge $erds_plus $gap_file

if test  -z ${cnv_dir:?"missing dir root"}
        then
                exit -1
        fi

erds_dir=$cnv_dir/erds
cnvn_dir=$cnv_dir/cnvn
merge_dir=$cnv_dir/erds+
alt_erds="alt.erds.txt"
alt_cnvn="alt.cnvn.txt"
cnv_tag=$cnv_dir/ALL_CNV

echo "Set-up.."
if  [[ -d $erds_dir/original ]] ; then
    echo "Found erds/original, creating erds/formatted"
    if  [[ ! -d $erds_dir/formatted ]] ; then
        mkdir $erds_dir/formatted $erds_dir/formatted/temp $erds_dir/merged/
    else
	rm $erds_dir/formatted/*.vcf $erds_dir/formatted/temp/* $erds_dir/merged/*
    fi
else
    echo "erds/original missing exit.."
    exit
fi

if  [[ -d $cnvn_dir/original ]] ; then
    echo "Found cnvn/original, creating cnvn/formatted_filtered"
    if  [[ -d $cnvn_dir/formatted_filtered ]] ; then
        if  [[ ! -d $cnvn_dir/formatted_filtered/formatted ]] ; then
	    mkdir $cnvn_dir/formatted_filtered/formatted $cnvn_dir/formatted_filtered/filtered $cnvn_dir/merged/
	else
            rm $cnvn_dir/formatted_filtered/formatted/* $cnvn_dir/formatted_filtered/filtered/* $cnvn_dir/merged/* 
        fi

    else
        mkdir $cnvn_dir/formatted_filtered $cnvn_dir/formatted_filtered/formatted $cnvn_dir/formatted_filtered/filtered $cnvn_dir/merged/
    fi
else
    echo "cnvn/original missing exit.."
    exit
fi

if  [[ -d $merge_dir ]] ; then
    echo "Found erds+, deleting files"
    rm -rf $merge_dir/* 
else
    mkdir $merge_dir
fi

echo "Formatting, filtering and merging cnvn output"
for f in `ls $cnvn_dir/original/`; do $python_path "$cnvn_format" $cnvn_dir/original/$f $cnvn_dir/formatted_filtered/$f; done
mv $cnvn_dir/formatted_filtered/*.calls.filtered.txt $cnvn_dir/formatted_filtered/filtered/
mv $cnvn_dir/formatted_filtered/*.txt $cnvn_dir/formatted_filtered/formatted/
for f in `ls $cnvn_dir/formatted_filtered/filtered/`; do g=`echo $f| sed 's/.calls.filtered.txt//g'`; echo $f"	"$g; done > $alt_cnvn
$python_path "$cnvn_merge" -i $cnvn_dir/formatted_filtered/formatted/ -a $alt_cnvn -o $cnvn_dir/merged/ -g "$gap_file" 

echo "Formatting and merging erds output"
for f in `ls $erds_dir/original/`; do $python_path "$erds_format" $erds_dir/original/$f $erds_dir/formatted/$f; done
mv $erds_dir/formatted/*.temp $erds_dir/formatted/temp
for f in `ls $erds_dir/formatted/ | grep -v temp`; do g=`echo $f| sed 's/.erds.vcf//g'`; echo $f"	"$g; done > $alt_erds
$python_path "$erds_merge" -i $erds_dir/formatted/ -a $alt_erds -o $erds_dir/merged/ -g "$gap_file"

echo "genearting sample ids.."
ls $erds_dir/original | sed 's/.erds.vcf//g' > ids.txt

echo "merging erds and cnvn calls"
for f in `cat ids.txt `; do $python_path "$erds_plus" -i $erds_dir/merged/$f".erds.vcf.cluster.txt" -a $cnvn_dir/merged/$f".calls.txt.cluster.txt" -o $merge_dir/$f.erds+.txt -s $f -c 0 -p reciprocal; done

echo "check log files.."

