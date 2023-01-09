#!/bin/bash

source ~/miniconda3/etc/profile.d/conda.sh

BAM=`file *.bam | cut -d: -f1`
BAI=`file *.bam.bai | cut -d: -f1`
outdir=`pwd`
input_bam="${outdir}/${BAM}"
bai_path="${outdir}/${BAI}"
id_name=`echo ${BAM} | cut -d. -f1`
TMP='False'

runIchorCNA='/mnt/hdd/public/share_resource/OS_WES/ichorCNA/scripts/runIchorCNA.R'

GC_19='/mnt/hdd/public/share_resource/OS_WES/ichorCNA/inst/extdata/gc_hg19_1000kb.wig' 
MAP_19='/mnt/hdd/public/share_resource/OS_WES/ichorCNA/inst/extdata/map_hg19_1000kb.wig'
GC_38='/mnt/hdd/public/share_resource/OS_WES/ichorCNA/inst/extdata/gc_hg38_1000kb.wig' 
MAP_38='/mnt/hdd/public/share_resource/OS_WES/ichorCNA/inst/extdata/map_hg38_1000kb.wig'
EXONs_19='/mnt/hdd/public/share_resource/OS_WES/ichorCNA/inst/extdata/Exon_regions_hg19.v0.bed'
EXONs_38='/mnt/hdd/public/share_resource/OS_WES/ichorCNA/inst/extdata/Exon_regions_hg38.v1.bed'


Help() {
    echo "usage: ./Auto-IchorCNA.sh [-i|n|s|r|c|o|t|h]"
    echo "optional argruments:"
    echo "    -i              Input Tumor file path (BAM file) *Default [all BAM in the current directory]"
    echo "                    If you fill out the BAM path, you should have the .bam.bai file in the same"
    echo "                    directory. If you choose the default value, this script will check for the"
    echo "                    existence of .bam.bai and generate it if it does not exist."
    echo "    -b              Input normal file path (BAM file) *In case you want to compare with normal."
    echo "                    This file should NOT where in same directory with Tumor BAM."
    echo "    -n              Id name of patient *Default [same as BAM file name]"
    echo "    -s              The type of next generation sequencing methods ('WGS' or 'WES') *Required"
    echo "    -r              Reference genome type ('hg19' or 'hg38') *Required"
    echo "    -c              chromosome name ('num' or 'chr') *Default [chr]"
    echo "    -o              Output directory path *Default [current directory]"
    echo "    -t              Temporary output (True or False) *Default [False]"
    echo "    -h              Show this help message and exit"
    echo
    echo "author's email:"
    echo "    songphon_sutthittha@cmu.ac.th"
    echo
    echo "** Please contact us if you have any questions or problems with this script."
    echo "------------------------------------------------------------------------------------------"
}

while getopts ":hi:n:r:o:t:s:c:b:" option; do
   case $option in
      h) # display Help
         Help
         exit;;
      i) # default [Current directory]
         input_bam=${OPTARG};;
      b) # default [Current directory]
         normal_bam=${OPTARG};;
      n) # default [BAM file name]
         id_name=${OPTARG};;
      s) # required ["WGS" or "WES"]
         sequencing_type=${OPTARG};;
      r) # required ["hg19" or "hg38"]
         reference_genome=${OPTARG};;
      c) # required ["hg19" or "hg38"]
         chromosome=${OPTARG};;
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
    TMP_path="${outdir}/TMP"
    mkdir ${TMP_path}
fi

#Check BAM.bai and Generate if it NOT exist
if ! [[ -e "${bai_path}" ]]; then
    if [[ "${TMP}" = 'True' ]]; then
        echo "BAI file is NOT exist" >> ${TMP_path}/01-Check_BAI.log
        echo "Start to generate BAI file" >> ${TMP_path}/01-Check_BAI.log
    fi
    module load samtools
    samtools index ${input_bam} ${input_bam}.bai
    if [[ "${TMP}" = 'True' ]]; then
        echo "BAI file Generation already finished" >> ${TMP_path}/01-Check_BAI.log
    fi
else
    if [[ "${TMP}" = 'True' ]]; then
        echo "BAI file is exist" >> ${TMP_path}/01-Check_BAI.log
    fi

fi

#Check BAM.bai of NORMAL and Generate if it NOT exist
if ! [[ -z "${normal_bam}" ]];then
    if [[ -e "${normal_bam}" ]];then
        normal_bai="${normal_bam}.bai"
        if ! [[ -e "${normal_bai}" ]]; then
            if [[ "${TMP}" = 'True' ]]; then
                echo "BAI file is NOT exist" >> ${TMP_path}/01-Check_BAI.log
                echo "Start to generate BAI file" >> ${TMP_path}/01-Check_BAI.log
            fi
            module load samtools
            samtools index ${normal_bam} ${normal_bam}.bai
            if [[ "${TMP}" = 'True' ]]; then
                echo "BAI file Generation already finished" >> ${TMP_path}/01-Check_BAI.log
            fi
        else
            if [[ "${TMP}" = 'True' ]]; then
                echo "BAI file is exist" >> ${TMP_path}/01-Check_BAI.log
            fi

        fi
    else
        echo -e "Error: normal BAM: ${normal_bam} is NOT exist"
    fi
fi


#chromosome type
if [[ "${chromosome}" = "chr" ]];then
    chromosome='chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY'
elif [[ "${chromosome}" = "num" ]];then
    chromosome='1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y'
fi

#Read Count bin Using HMMcopy with CONDA
ReadCounter() {
    if [[ "${TMP}" = 'True' ]]; then
        echo "Start to count read bins" >> ${TMP_path}/02-ReadCounter.log
    fi
    conda activate HMMcopy
    readCounter --window 1000000 \
    --quality 20 \
    --chromosome ${chromosome} \
    ${input_bam} > ${id_name}.wig 
    if [[ "${TMP}" = 'True' ]]; then
        echo "ReadCounter process already finished" >> ${TMP_path}/02-ReadCounter.log
    fi
    conda deactivate

    if ! [[ -z "${normal_bam}" ]];then
        readCounter --window 1000000 \
        --quality 20 \
        --chromosome ${chromosome} \
        ${normal_bam} > ${normal_bam}.wig
    fi
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
    conda deactivate
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
    conda deactivate
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
    conda deactivate
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
    conda deactivate
}

IchorCNA_WES_hg38_norm() {
    conda activate R_ichorCNA
    Rscript ${runIchorCNA} --id ${id_name} \
    --WIG ${id_name}.wig --ploidy "c(2,3)" \
    --NORMWIG ${normal_wig}.wig\
    --gcWig ${GC_38} \
    --mapWig ${MAP_38} \
    --exons.bed ${EXONs_38} \
    --includeHOMD False --chrs "c(1:22)" --chrTrain "c(1:22)" --genomeBuild "hg38" --genomeStyle "UCSC" \
    --estimateNormal True --estimatePloidy True --estimateScPrevalence True \
    --scStates "c(1,3)" --txnE 0.9999 --txnStrength 10000 --outDir ${outdir}
    conda deactivate
}


#Check BAM exist
if [[ -e "${input_bam}"  ]]; then
    echo "Input: ${input_bam}"
    if [[ "${sequencing_type}" = "WGS" ]];then
        if [[ "${reference_genome}" = "hg19" ]];then
            ReadCounter 
            if [[ "${TMP}" = 'True' ]]; then
                echo "Start to run IchorCNA program" >> ${TMP_path}/03-IchorCNA_Running.log
                IchorCNA_WGS_hg19 >> ${TMP_path}/03-IchorCNA_Running.log
            else
                IchorCNA_WGS_hg19

            fi
            echo
            echo "**All Process already finished"
            echo "------------------------------------------------------------------------------------------"
        elif [[ "${reference_genome}" = "hg38" ]];then
            ReadCounter 
            if [[ "${TMP}" = 'True' ]]; then
                echo "Start to run IchorCNA program" >> ${TMP_path}/03-IchorCNA_Running.log
                IchorCNA_WGS_hg38 >> ${TMP_path}/03-IchorCNA_Running.log
            else
                IchorCNA_WGS_hg38

            fi
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
            ReadCounter 
            if [[ "${TMP}" = 'True' ]]; then
                echo "Start to run IchorCNA program" >> ${TMP_path}/03-IchorCNA_Running.log
                IchorCNA_WES_hg19 >> ${TMP_path}/03-IchorCNA_Running.log
            else
                IchorCNA_WES_hg19 

            fi
            echo
            echo "**All Process already finished"
            echo "------------------------------------------------------------------------------------------"
        elif [[ "${reference_genome}" = "hg38" ]];then
            ReadCounter 
            if [[ "${TMP}" = 'True' ]]; then
                echo "Start to run IchorCNA program" >> ${TMP_path}/03-IchorCNA_Running.log
                IchorCNA_WES_hg38 >> ${TMP_path}/03-IchorCNA_Running.log
            else
                IchorCNA_WES_hg38 

            fi
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
    elif ! [[ -z "${normal_bam}" ]];then
        ReadCounter 
            if [[ "${TMP}" = 'True' ]]; then
                echo "Start to run IchorCNA program" >> ${TMP_path}/03-IchorCNA_Running.log
                IchorCNA_WES_hg38_norm >> ${TMP_path}/03-IchorCNA_Running.log
            else
                IchorCNA_WES_hg38_norm 

            fi
            echo
            echo "**All Process already finished"
            echo "------------------------------------------------------------------------------------------"
    else
        echo "Error: Unrecognized: ${sequencing_type} in arguments of sequencing method"
        exit
    fi
else
    echo "Error: ${input_bam} is NOT exists."
    exit
fi

