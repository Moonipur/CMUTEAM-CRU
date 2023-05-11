#!/bin/bash

module load gatk
module load samtools
module load picard


nextflow run GATK_preprocess.nf -with-dag preprocess_flowchart.png -with-trace -with-report report.html -resume
