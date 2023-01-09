#!/bin/bash

BAM=`file *.bam | cut -d: -f1`
outdir=`pwd`
input_bam="${outdir}/${BAM}"
id_name=`echo ${BAM} | cut -d. -f1`
TMP='False'
runIchorCNA='/mnt/hdd/public/share_resource/OS_WES/ichorCNA/scripts/runIchorCNA.R'

GC_19='/mnt/hdd/public/share_resource/OS_WES/ichorCNA/inst/extdata/gc_hg19_1000kb.wig' 
MAP_19='/mnt/hdd/public/share_resource/OS_WES/ichorCNA/inst/extdata/map_hg19_1000kb.wig'
GC_38='/mnt/hdd/public/share_resource/OS_WES/ichorCNA/inst/extdata/gc_hg38_1000kb.wig' 
MAP_38='/mnt/hdd/public/share_resource/OS_WES/ichorCNA/inst/extdata/map_hg38_1000kb.wig'
EXONs_19='/mnt/hdd/public/share_resource/OS_WES/ichorCNA/inst/extdata/Exon_regions_hg19.v0.bed'
EXONs_38='/mnt/hdd/public/share_resource/OS_WES/ichorCNA/inst/extdata/Exon_regions_hg38.v1.bed'

eval "$(conda shell.bash hook)"

Help() {
    echo "usage: ./Auto-IchorCNA.sh [-i|n|s|r|o|t|h]"
    echo "optional argruments:"
    echo "    -i              input file path (BAM file) *Default [all BAM in the current directory]"
    echo "    -n              id name of patient *Default [same as BAM file name]"
    echo "    -s              the type of next generation sequencing methods ("WGS" or "WES") *Required"
    echo "    -r              reference genome type ("hg19" or "hg38") *Required"
    echo "    -o              output directory path *Default [current directory]"
    echo "    -t              temporary output (True or False) *Default [False]"
    echo "    -h              show this help message and exit"
    echo
    echo "author's email:"
    echo "    songphon_sutthittha@cmu.ac.th"
    echo
    echo "** Please contact us if you have any questions or problems with this script."
    echo "------------------------------------------------------------------------------------------"
}

while getopts ":hi:n:r:o:t:s:" option; do
   case $option in
      h) # display Help
         Help
         exit;;
      i) # default [Current directory]
         input_bam=${OPTARG};;
      n) # default [BAM file name]
         id_name=${OPTARG};;
      s) # required ["WGS" or "WES"]
         sequencing_type=${OPTARG};;
      r) # required ["hg19" or "hg38"]
         reference_genome=${OPTARG};;
      o) # default [Current directory]
         outdir=${OPTARG};;
      t) # default [False]
         TMP=${OPTARG};;
     \?) # Invalid option
         echo "Error: Unrecognized arguments"
         exit;;
   esac
done

#Check Output directory exist
if [[ -e "${outdir}" ]]; then
    # Remove / at the last character of path
    if [[ "${outdir: -1}" = "/" ]]; then
        outdir=`echo ${outdir} | sed 's/.$//'`
        echo "Output: ${outdir}"
    else
        echo "Output: ${outdir}"
    fi
else
    echo "Error: ${outdir} is NOT exists."
    exit
fi


#Create temporary directory
if [[ "${TMP}" = 'True' ]]; then
    mkdir ${outdir}/TMP
fi

#Read Count bin Using HMMcopy with CONDA
ReadCounter() {
    conda activate HMMcopy
    readCounter --window 1000000 \
    --quality 20 \
    --chromosome "1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y" \
    ${input_bam} > ${id_name}.wig
    conda deactivate
}

#Run ichorCNA
IchorCNA_WGS_hg19() {
    conda activate R_ichorCNA
    Rscript ${runIchorCNA} --id ${id_name} \
    --WIG ${id_name}.wig --ploidy "c(2,3)" \
    --gcWig ${GC_19} \
    --mapWig ${MAP_19} \
    --includeHOMD False --chrs "c(1:22)" --chrTrain "c(1:22)" --genomeBuild "hg19" --genomeStyle "UCSC" \
    --estimateNormal True --estimatePloidy True --estimateScPrevalence True \
    --scStates "c(1,3)" --txnE 0.9999 --txnStrength 10000 --outDir ${outdir}
}

IchorCNA_WGS_hg38() {
    conda activate R_ichorCNA
    Rscript ${runIchorCNA} --id ${id_name} \
    --WIG ${id_name}.wig --ploidy "c(2,3)" \
    --gcWig ${GC_38} \
    --mapWig ${MAP_38} \
    --includeHOMD False --chrs "c(1:22)" --chrTrain "c(1:22)" --genomeBuild "hg38" --genomeStyle "UCSC" \
    --estimateNormal True --estimatePloidy True --estimateScPrevalence True \
    --scStates "c(1,3)" --txnE 0.9999 --txnStrength 10000 --outDir ${outdir}
}

IchorCNA_WES_hg19() {
    conda activate R_ichorCNA
    Rscript ${runIchorCNA} --id ${id_name} \
    --WIG ${id_name}.wig --ploidy "c(2,3)" \
    --gcWig ${GC_19} \
    --mapWig ${MAP_19} \
    --exons.bed ${EXONs_19} \
    --includeHOMD False --chrs "c(1:22)" --chrTrain "c(1:22)" --genomeBuild "hg19" --genomeStyle "UCSC" \
    --estimateNormal True --estimatePloidy True --estimateScPrevalence True \
    --scStates "c(1,3)" --txnE 0.9999 --txnStrength 10000 --outDir ${outdir}
}

IchorCNA_WES_hg38() {
    conda activate R_ichorCNA
    Rscript ${runIchorCNA} --id ${id_name} \
    --WIG ${id_name}.wig --ploidy "c(2,3)" \
    --gcWig ${GC_38} \
    --mapWig ${MAP_38} \
    --exons.bed ${EXONs_38} \
    --includeHOMD False --chrs "c(1:22)" --chrTrain "c(1:22)" --genomeBuild "hg38" --genomeStyle "UCSC" \
    --estimateNormal True --estimatePloidy True --estimateScPrevalence True \
    --scStates "c(1,3)" --txnE 0.9999 --txnStrength 10000 --outDir ${outdir}
}


#Check BAM exist
if [[ -e "${input_bam}"  ]]; then
    echo "Input: ${input_bam}"
    if [[ "${sequencing_type}" = "WGS" ]];then
        if [[ "${reference_genome}" = "hg19" ]];then
            IchorCNA_WGS_hg19
            echo
            echo "**All Process already finished"
            echo "------------------------------------------------------------------------------------------"
        elif [[ "${reference_genome}" = "hg38" ]];then
            IchorCNA_WGS_hg38
            echo
            echo "**All Process already finished"
            echo "------------------------------------------------------------------------------------------"
        elif [[ -z "${reference_genome}" ]];then
            echo "Error: Reference human genome type is required. **Please fill out"
            exit
        else
            echo "Error: Unrecognized: ${reference_genome} in arguments of reference genome"
            exit
        fi
    elif [[ "${sequencing_type}" = "WES" ]];then
        if [[ "${reference_genome}" = "hg19" ]];then
            IchorCNA_WES_hg19
            echo
            echo "**All Process already finished"
            echo "------------------------------------------------------------------------------------------"
        elif [[ "${reference_genome}" = "hg38" ]];then
            IchorCNA_WES_hg38
            echo
            echo "**All Process already finished"
            echo "------------------------------------------------------------------------------------------"
        elif [[ -z "${reference_genome}" ]];then
            echo "Error: Reference human genome type is required. **Please fill out"
            exit
        else
            echo "Error: Unrecognized: ${reference_genome} in arguments of reference genome"
            exit
        fi
    elif [[ -z "${sequencing_type}" ]];then
        echo "Error: Sequencing method is required. **Please fill out"
        exit
    else
        echo "Error: Unrecognized: ${sequencing_type} in arguments of sequencing method"
        exit
    fi
else
    echo "Error: ${input_bam} is NOT exists."
    exit
fi

