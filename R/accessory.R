#' @export 

GRanges <- function(seqnames, ranges, gc, mappability) {
    GenomicRanges::GRanges(
        seqnames = seqnames, 
        ranges = ranges, 
        gc = gc,
        mappability = mappability
    )
}

#' @export 
IRanges <- function(start, end) {
    IRanges::IRanges(start = start, end = end)
}

#' @export 
accessory.file_exists <- function(allfile_path) {
    all_exists = TRUE
    for(file in allfile_path) {
        if(!file.exists(file)) {
            paste(file,"does NOT exists", sep=" ")
            all_exists = FALSE
        } else {
            paste(file,"exists", sep=" ")
        }
    }
    return(all_exists)
}

#' @export 
check_bai_path <- function(bam_path) {
    bai_path = paste0(bam_path,".bai")
    if(file.exists(bai_path)) {
        return(TRUE)
    } else {
        return(FALSE)
    }
}

#' @export 
check_gzip <- function(target_bed){
  if(summary( file(target_bed) )$class == "gzfile")
    TRUE
  else FALSE
}

#' @export 
check_Chr_UCSU_format <- function(bam_path) {
    all_chrs = Rsamtools::idxstatsBam(bam_path)$seqnames
    return(any(all_chrs %in% c("chr1")))
}

#' @export 
convert_GRch2UCSC <- function(which) {
    GenomeInfoDb::seqlevels(which) <- sub('chrM', 'M', GenomeInfoDb::seqlevels(which))
    GenomeInfoDb::seqlevels(which) <- gsub('^(.*)$', 'chr\\1', GenomeInfoDb::seqlevels(which))
    return(which)
}

#' @export 
convert_UCSC2Chr <- function(bam) {
  names(bam) = gsub(
    "^chr","", names(bam)
  )

  return(bam)
}

#' @export 
SAVE_toCSV <- function(data, path, index=FALSE) {
    write.csv(data, path, row.names=index)
}