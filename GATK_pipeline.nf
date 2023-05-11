/* 
 * Nextflow version 21.10.6
 * usage: nextflow run GATK_preprocess.nf -with-dag preprocess_flowchart.png -with-trace -with-report report.html -resume
 * nextflow self-update
 * nextflow.enable.dsl=2
 */
 
nextflow.enable.dsl=1

params.reads = "/mnt/sas/public/share_resource/OS_WES/cfDNA/OS3_R{1,2}_trimmed.fastq.gz"
params.outdir = "OS03"
params.knownSites_dbsnp="/mnt/sas/ref/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz"
params.knownSites_indel="/mnt/sas/ref/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
params.ref="/mnt/sas/ref/hg38/v0/Homo_sapiens_assembly38.fasta"
params.memory=10
params.cpus=16


log.info """\
         C M U T E A M - N F   P I P E L I N E    
         =====================================
         pipeline name		: fastq preprocess based on GATK best practice
         pipeline version	: 0.1.1
         maintainer		: Apiwat Sangphukieo 
         contact		: apiwat.sang@cmu.ac.th
         edited by		: Songphon Sutthiithasakul
         developed under	: Nextflow version 21.10.6         
         reads			: ${params.reads}
         outdir			: ${params.outdir}
         knownSites_dbsnp	: ${params.knownSites_dbsnp}
         knownSites_indel	: ${params.knownSites_indel}
         ref			: ${params.ref}
         memory			: ${params.memory}
         cpus			: ${params.cpus}
         """
         .stripIndent()

// check if the file exists
if( !file(params.knownSites_dbsnp).exists() ) exit 1, "Missing knownSites dbsnp file: Homo_sapiens_assembly38.dbsnp138.vcf.gz"
if( !file(params.knownSites_indel).exists() ) exit 1, "Missing knownSites indel file: Homo_sapiens_assembly38.known_indels.vcf.gz"
if( !file(params.ref).exists() ) exit 1, "Missing ref file: Homo_sapiens_assembly38.fasta"
if( !file(params.ref+".fai").exists() ) exit 1, "Missing index ref file: Homo_sapiens_assembly38.fasta.fai"
if( !file(params.ref.replace(".fasta",".dict")).exists() ) exit 1, "Missing index ref file: Homo_sapiens_assembly38.dict"

ref_fasta=file(params.ref)
ref_fasta_fai = file("${params.ref}.fai")
ref_dict = file(params.ref.replace(".fasta",".dict") )

known_indel = file(params.knownSites_indel)
known_indel_index = file("${params.knownSites_indel}.tbi")
known_dbSNP = file(params.knownSites_dbsnp)
known_dbSNP_index = file("${params.knownSites_dbsnp}.tbi")


Channel 
    .fromFilePairs( params.reads, checkIfExists: true )
    .into { read_pairs_ch; read_pairs2_ch} 


process fastqc {
    tag "FASTQC on $pair_id"
    publishDir "$params.outdir/$pair_id/01_FASTQC", mode:'copy'

    input:
    set val(pair_id), file(reads) from read_pairs_ch

    output:
    set val(pair_id), file("fastqc_${pair_id}_logs") into fastqc_ch

    script:
    """
    mkdir fastqc_${pair_id}_logs
    fastqc -o fastqc_${pair_id}_logs -f fastq -q ${reads} &
    """  
}  
 

process BWA_MEM {
    cpus "${params.cpus}"
    memory "${params.memory}"
    tag "BWA MEM on $pair_id"
    publishDir "$params.outdir/$pair_id/03_BWA_READ_ALIGNMENT", mode:'copy'   

    input:
    set pair_id, file(reads) from read_pairs2_ch

    
    path(ref_dict)
    path(ref_fasta)
    path(ref_fasta_fai)

    output:
      set val(pair_id) ,file("${pair_id}.bam") into bam_ch, bam_ch2
      file("software_versions.${task.process}.txt")

    """
    bwa mem -R '@RG\\tID:${pair_id}\\tPL:ILLUMINA\\tSM:${pair_id}' -t 10 ${params.ref} ${reads[0]} ${reads[1]} | samtools view -@${task.cpus} -bS - > ${pair_id}.bam
    echo "bwa=0.7.17"  >> software_versions.${task.process}.txt
    samtools --version >> software_versions.${task.process}.txt
    """
}

process FRAGSTAT {
    tag "FRAGSTAT on $pair_id"
    publishDir "$params.outdir/$pair_id/03_BWA_READ_ALIGNMENT", mode:'copy'  

    input:
    set pair_id, file(bam) from bam_ch2

    output:
      file("${pair_id}_alignment_metrics.txt")

    """
    samtools flagstat ${pair_id}.bam > ${pair_id}_alignment_metrics.txt
    """
}

process PICARD_SORTSAM {
    memory "${params.memory}"
    tag "PICARD_SORTSAM on $pair_id"
    publishDir "$params.outdir/$pair_id/04_PICARD_SORTSAM", mode:'copy'  

    input:
    set pair_id, file(bam) from bam_ch

    output:
      set val(pair_id), file("${pair_id}_sorted.bam") into sort_bam_ch
      file("${pair_id}_sorted.bai") into sort_bai_ch

    """
    java -Xmx${params.memory}G -jar $PICARD_JAR SortSam \
	I=${bam} \
	O=${pair_id}_sorted.bam \
	SORT_ORDER=coordinate \
	TMP_DIR=${baseDir}/TMP \
	CREATE_INDEX=true
   
    """
}

process MARK_DUP {
    memory "${params.memory}"
    tag "MARK_DUP on $pair_id"
    publishDir "$params.outdir/$pair_id/05_MARK_DUPLICATE", mode:'copy'  

    input:
    set pair_id, file(bam) from sort_bam_ch
    file(bai) from sort_bai_ch

    output:
      set val(pair_id), file("${pair_id}_sorted_dedup.bam") into dedup_bam_ch, dedup_bam_ch2
      file("${pair_id}_mark_dup_metrics.txt")

    """
    java -jar -Xmx${params.memory}G $PICARD_JAR MarkDuplicates \
	INPUT=${bam} \
	OUTPUT=${pair_id}_sorted_dedup.bam \
	METRICS_FILE=${pair_id}_mark_dup_metrics.txt \
	TMP_DIR=${baseDir}/TMP
   
    """
}

process BaseRecalibrator {
    memory "${params.memory}"
    tag "BaseRecalibrator on $pair_id"
    publishDir "$params.outdir/$pair_id/06_BaseRecalibrator", mode:'copy'  

    input:
    set pair_id, file(dedup_bam) from dedup_bam_ch
    path(ref_fasta)
    path(ref_fasta_fai)
    path(ref_dict)
    path(known_indel)
    path(known_dbSNP)
    path(known_dbSNP_index)
    path(known_indel_index)
    
    output:
      file("${pair_id}_recal.table") into BaseRecal_ch
      

    """
    gatk --java-options "-Xms8G -Xmx${params.memory}G" BaseRecalibrator \
	-R ${ref_fasta} \
	-I ${dedup_bam} \
	-known-sites ${known_dbSNP} \
	-known-sites ${known_indel} \
	-O ${pair_id}_recal.table
    """
}


process ApplyBQSR {
    memory "${params.memory}"
    tag "ApplyBQSR on $pair_id"
    publishDir "$params.outdir/$pair_id/07_ApplyBQSR", mode:'copy'  

    input:
	set pair_id, file(dedup_bam) from dedup_bam_ch2
	file(table) from BaseRecal_ch
	path(ref_dict)
	path(ref_fasta)
	path(ref_fasta_fai)
    
    output:
      file("${pair_id}_recal.bam") into BQSR_bam_ch
      file("${pair_id}_recal.bai") into BQSR_bai_ch

    """
    gatk --java-options "-Xms8G -Xmx${params.memory}G" ApplyBQSR \
	-R ${ref_fasta} \
	-I ${dedup_bam} \
	-bqsr ${table} \
	-O ${pair_id}_recal.bam
    """
}

workflow.onComplete { 
	log.info ( workflow.success ? "\nDone! Open the following QC report in your browser --> $params.outdir/multiqc_report.html\n" : "Oops .. something went wrong" )
}
