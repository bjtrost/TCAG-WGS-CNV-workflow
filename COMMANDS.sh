# This document contains examples of commands used to run BWA and the CNV-detection algorithms.
# Of course, details will differ depending on how your system is set up, so these commands
# are intended only as a guide.

### BWA
bwa mem -M genome.fa read1.fastq read2.fastq > output.sam

### GATK
# Followed standard best practices recommendations

### Canvas
$PREFIX/mytools/mono/mono-4.0.2/bin/mono $PREFIX/mytools/canvas/canvas-1.3.5_x64/Canvas.exe Germline-WGS \\
        -b ${bamfile} \\
        --b-allele-vcf ${vcffile} \\
        -o ${work_dir} \\
        -r ${CANVASPREFIX}/kmer.fa \\
        -g ${CANVASPREFIX}/WholeGenomeFasta \\
        -f ${CANVASPREFIX}/filter13.bed \\
        -n ${sample}

# cn.MOPS (R code)
bamDataRanges <- getReadCountsFromBAM(BAMFiles, refSeqName=chrom, mode='paired')
resCNMOPS <- cn.mops(bamDataRanges)
resCNMOPS <- calcIntegerCopyNumbers(resCNMOPS)
cnvr(resCNMOPS)

### CNVnator
cnvnator -root ${rootfile} -genome $hg  -unique -tree ${bamfile[*]}
cnvnator -root ${rootfile} -genome $hg  -d `dirname $genome` -his $binsize
cnvnator -root ${rootfile} ${chromosome} -stat $binsize
cnvnator -root ${rootfile} ${chromosome} -partition $binsize
cnvnator -root ${rootfile} ${chromosome} -call $binsize

### ERDS
perl $ERDSPATH/erds_pipeline.pl -o ${work_dir}/${samplename} -b ${used_bamfile} -v ${used_vcffile} -r $genome

### GenomeSTRiP
# The steps below are modifications of their sample script, tailored for our system (PBS job system)
# You may need to follow different steps if you use a different system.
# Since GenomeSTRiP does not support the PBS job system, we installed Drmaa, which simulates the SGE job system supported by GenomeSTRiP (-jobRunner Drmaa)

# Step1
java -cp ${classpath} ${mx} \
        org.broadinstitute.gatk.queue.QCommandLine \
        -S ${SV_DIR}/qscript/SVPreprocess.q \
        -S ${SV_DIR}/qscript/SVQScript.q \
        -gatk ${SV_DIR}/lib/gatk/GenomeAnalysisTK.jar \
        --disableJobReport \
        -cp ${classpath} \
        -configFile ${configfile} \
        -tempDir ${SV_TMPDIR} \
        -R ${used_genome} \
        -genomeMaskFile ${g_masked_genome} \
        -copyNumberMaskFile ${cnv_masked_genome} \
        -ploidyMapFile ${ploidy_genome} \
        -genderMapFile ${genderfile} \
        -md ${preprocessDir}/metadata \
        -runDirectory ${preprocessDir} \
        -disableGATKTraversal \
        -useMultiStep \
        -reduceInsertSizeDistributions false \
        -computeGCProfiles true \
        -computeReadCounts true \
        -jobLogDir ${preprocessDir}/logs \
        ${bamfiles} \
        -bamFilesAreDisjoint true\
        -jobRunner Drmaa \
        -jobProject ${projname} \
        -jobNative " -l vmem=32gb" \
        -debug true \
        -run || exit 1

# Step 2
java -cp ${classpath} ${mx} \
    org.broadinstitute.gatk.queue.QCommandLine \
    -S ${SV_DIR}/qscript/SVDiscovery.q \
    -S ${SV_DIR}/qscript/SVQScript.q \
    -gatk ${SV_DIR}/lib/gatk/GenomeAnalysisTK.jar \
    --disableJobReport \
    -cp ${classpath} \
    -configFile ${configfile} \
    -tempDir ${SV_TMPDIR} \
    -R ${used_genome}\
    -genomeMaskFile ${cnv_masked_genome} \
    -ploidyMapFile ${ploidy_genome} \
    -genderMapFile ${genderfile} \
    -md ${preprocessDir}/metadata \
    -runDirectory ${discoveryDir} \
    -disableGATKTraversal \
    -jobLogDir ${discoveryDir}/logs \
    -minimumSize 100 \
    -maximumSize 1000000 \
    -suppressVCFCommandLines \
    ${bamfiles} \
    -O ${sites} \
    -jobRunner Drmaa \
    -jobProject ${projname} \
    -jobNative " -l vmem=8gb " \
    -P select.validateReadPairs:false \
    -run || exit 1

# Step 3
java -cp ${classpath} ${mx} \
    org.broadinstitute.gatk.queue.QCommandLine \
    -S ${SV_DIR}/qscript/SVGenotyper.q \
    -S ${SV_DIR}/qscript/SVQScript.q \
    -gatk ${SV_DIR}/lib/gatk/GenomeAnalysisTK.jar \
    --disableJobReport \
    -cp ${classpath} \
    -configFile ${configfile} \
    -tempDir ${SV_TMPDIR} \
    -R ${used_genome} \
    -genomeMaskFile ${cnv_masked_genome} \
    -ploidyMapFile ${ploidy_genome}  \
    -genderMapFile ${genderfile}     \
    -runDirectory ${genotypingDir}   \
    -md ${preprocessDir}/metadata    \
    -disableGATKTraversal \
    -jobLogDir ${preprocessDir}/logs \
    ${bamfiles}   \
    -vcf ${sites} \
    -O ${genotypes}  \
    -jobRunner Drmaa \
    -jobProject ${projname}    \
    -jobNative " -l vmem=8gb " \
    -P select.validateReadPairs:false \
    -debug true \
    -run


### RDXplorer
winSize=100
baseCopy=2
filter=10
sumWithZero=True
debug=True
delete=False

$RDXPLORER/rdxplorer.py \
       ${path2bam}        \
       ${reference}       \
       ${workdir}         \
       ${chromOfInterest} \
       ${gender}          \
       ${hg}              \
       ${winSize}         \
       ${baseCopy}        \
       ${filter}          \
       ${sumWithZero}     \
       ${debug}           \
       ${delete}
