#!/bin/bash

module load samtools

cur_dir=`pwd`
BAM=$1
count_slash=`echo ${BAM} | grep -o / | wc -l`
slash=$((${count_slash} + 1))
new_BAM=`echo ${BAM} | cut -d/ -f ${slash} | cut -d. -f 1`

#Check type of Sizeselection
if ![[ "$2" = "Both" ]]; then
    if [[ "$2" = "150" OR "$2" = "167" ]]; then
        size=$2
    else
        echo "Error: our program have NOT cut-off at size ($2), please check again."
        exit
    fi
else
    size="Both"
fi

#Check similar directory
if [[ $3 = "${cur_dir}/"]]; then
    out_dir="${cur_dir}/"
else
    out_dir=$3
fi

new_outdir="${out_dir}SizeSelection_out"

mkdir ${new_outdir}

#Read size selection
Sort150() {
    samtools view -h ${BAM} | \
    awk 'substr($0,1,1)=="@" || ($9 >= -150 && $9 <= 150 && $9 != 0)' | \
    samtools view -b > ${new_outdir}/${new_BAM}_sort150.bam
}

Sort167() {
    samtools view -h ${BAM} | \
    awk 'substr($0,1,1)=="@" || ($9 >= -167 && $9 <= 167 && $9 != 0)' | \
    samtools view -b > ${new_outdir}/${new_BAM}_sort167.bam
}

#Have limitation (can't screen read larger than 1000 bp)
More150() {
    samtools view -h ${BAM} | \
    awk 'substr($0,1,1)=="@" || ($9 < -150 && $9 >= -1000) || ($9 > 150 && $9 <= 1000)' | \
    samtools view -b > ${new_outdir}/${new_BAM}_more150.bam
}

More167() {
    samtools view -h ${BAM} | \
    awk 'substr($0,1,1)=="@" || ($9 < -167 && $9 >= -1000) || ($9 > 167 && $9 <= 1000)' | \
    samtools view -b > ${new_outdir}/${new_BAM}_more167.bam
}

if [[ ${size} = 150 ]]; then
    Sort150()
    More150()
elif [[ ${size} = 167 ]]; then
    Sort167()
    More167()
elif [[ ${size} = "Both" ]]; then
    Sort150()
    Sort167()
    More150()
    More167()
fi
