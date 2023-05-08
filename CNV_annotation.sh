#!/bin/bash

source ~/miniconda3/etc/profile.d/conda.sh
PACKAGE="/mnt/hdd/public/share_resource/OS_WES"
cur_dir=`pwd`

Help() {
    echo "usage: CNV_annotation [-i|-o|-r|h]"
    echo "description: CNV_annotation is a software suite script that uses the ACE and ClassifyCNV packages"
    echo "             for annotating gains and losses of CNVs and predicting list of dosage-sensitive genes"
    echo "             that affect as pathogenic variant."
    echo "optional argruments:"
    echo "    -i       Input file path (BAM file)"
    echo "    -o       Output directory path (/path/of/output/directory/), Default is the currect directory"
    echo "    -r       Type of reference human genome (hg19/hg38)"
    echo "    -h       Show this help message and exit"
    echo
    echo "author's email:"
    echo "    songphon_sutthittha@cmu.ac.th"
    echo
    echo "** Please contact us if you have any questions or problems with this script."
    echo "------------------------------------------------------------------------------------------"
}

while getopts ":hi:o:r:" option; do
   case $option in
      h) # display Help
         Help
         exit;;
      i) # default [Current directory]
         name=${OPTARG};;
      o) # default [Current directory]
         outdir=${OPTARG};;
      r) # default [Current directory]
         ref=${OPTARG};;
     \?) # Invalid option
         echo "Error: Unrecognized arguments"
         exit;;
   esac
done

NAME=`echo ${name}| cut -d. -f1`

if [ ! -z ${outdir} ]
then
    cd ${outdir}
    OUTDIR=`pwd`
    cd ${cur_dir}
else
    OUTDIR=`pwd`
fi

#Step 1
echo -e ">> 1ST: The ACE scipt is generating and reading using the ACE software"
conda activate Rdevtools
if [ "${ref}" = "hg19" ]
then
    if [ ! -f ${OUTDIR}/ACE_SCRIPT_hg19.R ]
    then
        echo "#!/usr/bin/env Rscript" > ${OUTDIR}/ACE_SCRIPT_hg19.R
        echo '' >> ${OUTDIR}/ACE_SCRIPT_hg19.R
        echo 'library(QDNAseq)' >> ${OUTDIR}/ACE_SCRIPT_hg19.R
        echo '' >> ${OUTDIR}/ACE_SCRIPT_hg19.R
        echo 'library(ACE)' >> ${OUTDIR}/ACE_SCRIPT_hg19.R
        echo 'bins <- getBinAnnotations(binSize=100, genome="hg19")' >> ${OUTDIR}/ACE_SCRIPT_hg19.R
        echo 'readCounts <- binReadCounts(bins)' >> ${OUTDIR}/ACE_SCRIPT_hg19.R
        echo 'readCountsFiltered <- applyFilters(readCounts, residual=TRUE, blacklist=TRUE)' >> ${OUTDIR}/ACE_SCRIPT_hg19.R
        echo 'readCountsFiltered <- estimateCorrection(readCountsFiltered)' >> ${OUTDIR}/ACE_SCRIPT_hg19.R
        echo 'copyNumbers <- correctBins(readCountsFiltered)' >> ${OUTDIR}/ACE_SCRIPT_hg19.R
        echo 'copyNumbersNormalized <- normalizeBins(copyNumbers)' >> ${OUTDIR}/ACE_SCRIPT_hg19.R
        echo 'copyNumbersSmooth <- smoothOutlierBins(copyNumbersNormalized)' >> ${OUTDIR}/ACE_SCRIPT_hg19.R
        echo 'copyNumbersSegmented <- segmentBins(copyNumbersSmooth, transformFun="sqrt")' >> ${OUTDIR}/ACE_SCRIPT_hg19.R
        echo 'object <- copyNumbersSegmented' >> ${OUTDIR}/ACE_SCRIPT_hg19.R
        echo 'template <- objectsampletotemplate(object,index=1)' >> ${OUTDIR}/ACE_SCRIPT_hg19.R
        echo 'segmentdf <- getadjustedsegments(template, cellularity = 0.25)' >> ${OUTDIR}/ACE_SCRIPT_hg19.R
        echo 'sqmodel <- squaremodel(template, prows = 150, ptop = 3.3, pbottom = 1.8, ' >> ${OUTDIR}/ACE_SCRIPT_hg19.R
        echo '                       penalty = 0.5, penploidy = 0.5)' >> ${OUTDIR}/ACE_SCRIPT_hg19.R
        echo '' >> ${OUTDIR}/ACE_SCRIPT_hg19.R
        echo 'models <- data.frame(segmentdf)' >> ${OUTDIR}/ACE_SCRIPT_hg19.R
        echo -e "write.table(models, file.path('${OUTDIR}/${NAME}_segmentation.tsv'), quote = FALSE, " >> ${OUTDIR}/ACE_SCRIPT_hg19.R
        echo "            sep = '\t', na = '', row.names = FALSE)" >> ${OUTDIR}/ACE_SCRIPT_hg19.R
        echo "" >> ${OUTDIR}/ACE_SCRIPT_hg19.R
        Rscript ${OUTDIR}/ACE_SCRIPT_hg19.R
    else
        if [ ! -f ${OUTDIR}/${NAME}_segmentation.tsv ]
        then
            Rscript ${OUTDIR}/ACE_SCRIPT_hg19.R
        fi
    fi
elif [ "${ref}" = "hg38" ]
then
    if [ ! -f ${OUTDIR}/ACE_SCRIPT_hg38.R ]
    then
        echo "#!/usr/bin/env Rscript" > ${OUTDIR}/ACE_SCRIPT_hg38.R
        echo '' >> ${OUTDIR}/ACE_SCRIPT_hg38.R
        echo 'library(QDNAseq)' >> ${OUTDIR}/ACE_SCRIPT_hg38.R
        echo 'library(QDNAseq.hg38)' >> ${OUTDIR}/ACE_SCRIPT_hg38.R
        echo 'library(ACE)' >> ${OUTDIR}/ACE_SCRIPT_hg38.R
        echo 'bins <- getBinAnnotations(binSize=100, genome="hg38")' >> ${OUTDIR}/ACE_SCRIPT_hg38.R
        echo 'readCounts <- binReadCounts(bins)' >> ${OUTDIR}/ACE_SCRIPT_hg38.R
        echo 'readCountsFiltered <- applyFilters(readCounts, residual=TRUE, blacklist=TRUE)' >> ${OUTDIR}/ACE_SCRIPT_hg38.R
        echo 'readCountsFiltered <- estimateCorrection(readCountsFiltered)' >> ${OUTDIR}/ACE_SCRIPT_hg38.R
        echo 'copyNumbers <- correctBins(readCountsFiltered)' >> ${OUTDIR}/ACE_SCRIPT_hg38.R
        echo 'copyNumbersNormalized <- normalizeBins(copyNumbers)' >> ${OUTDIR}/ACE_SCRIPT_hg38.R
        echo 'copyNumbersSmooth <- smoothOutlierBins(copyNumbersNormalized)' >> ${OUTDIR}/ACE_SCRIPT_hg38.R
        echo 'copyNumbersSegmented <- segmentBins(copyNumbersSmooth, transformFun="sqrt")' >> ${OUTDIR}/ACE_SCRIPT_hg38.R
        echo 'object <- copyNumbersSegmented' >> ${OUTDIR}/ACE_SCRIPT_hg38.R
        echo 'template <- objectsampletotemplate(object,index=1)' >> ${OUTDIR}/ACE_SCRIPT_hg38.R
        echo 'segmentdf <- getadjustedsegments(template, cellularity = 0.25)' >> ${OUTDIR}/ACE_SCRIPT_hg38.R
        echo 'sqmodel <- squaremodel(template, prows = 150, ptop = 3.3, pbottom = 1.8, ' >> ${OUTDIR}/ACE_SCRIPT_hg38.R
        echo '                       penalty = 0.5, penploidy = 0.5)' >> ${OUTDIR}/ACE_SCRIPT_hg38.R
        echo '' >> ${OUTDIR}/ACE_SCRIPT_hg38.R
        echo 'models <- data.frame(segmentdf)' >> ${OUTDIR}/ACE_SCRIPT_hg38.R
        echo -e "write.table(models, file.path('${OUTDIR}/${NAME}_segmentation.tsv'), quote = FALSE, " >> ${OUTDIR}/ACE_SCRIPT_hg38.R
        echo "            sep = '\t', na = '', row.names = FALSE)" >> ${OUTDIR}/ACE_SCRIPT_hg38.R
        echo "" >> ${OUTDIR}/ACE_SCRIPT_hg38.R
        Rscript ${OUTDIR}/ACE_SCRIPT_hg38.R
    else
        if [ ! -f ${OUTDIR}/${NAME}_segmentation.tsv ]
        then
            Rscript ${OUTDIR}/ACE_SCRIPT_hg38.R
        fi
    fi
fi
conda deactivate

#Step 2
echo -e ">> 2ND: The BED file is converting from TSV file using Python module"
conda activate Bedtools
cut -f1,2,3,8 ${OUTDIR}/${NAME}_segmentation.tsv | sed -z 's/\t/,/g' > ${OUTDIR}/${NAME}_segmentation_cut.csv
python3 ${PACKAGE}/Code_session/Classify_Annotation/TSV2BED.py ${OUTDIR}/${NAME}_segmentation_cut.csv ${OUTDIR} ${NAME}_segmentation.csv
sed -z 's/,/\t/g' ${OUTDIR}/${NAME}_segmentation.csv > ${OUTDIR}/${NAME}_segmentation.bed
rm ${OUTDIR}/${NAME}_segmentation_cut.csv
rm ${OUTDIR}/${NAME}_segmentation.csv

#Step 3
echo -e ">> 3RD: The list of dosage-sensitive genes and pathogenic variants are predicting using ClassifyCNV package"
cd ${outdir}
python3 ${PACKAGE}/Code_session/ClassifyCNV/ClassifyCNV.py --infile ${NAME}_segmentation.bed --GenomeBuild ${ref}
cd ${cur_dir}
conda deactivate

echo -e ">> 4TH: All steps already finish and saved out at; ${OUTDIR}"