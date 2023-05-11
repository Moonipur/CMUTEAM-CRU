#' @export

import_bamfile <- function(
    bam_path, bin_size = 1000, ref_genome = "hg38",
    target_bed = NULL, min_mapQ = 20, blacklist = TRUE
) {
    ## check bamfile path exists
    if(!accessory.file_exists(bam_path)) {
        stop("Error: Your bamfile path does NOT exists. Please check bamfile path!")
    }

    ## check target bedfile exists
    if(!is.null(target_bed) & !accessory.file_exists(target_bed)) {
        stop("Error: Your target bedfile does NOT exists.")
    }

    ## check baifile exists and generating in case dose NOT exists
    if(!check_bai_path(bam_path)) {
        print("Your baifile is missing. Starting to generate baifile.")
        Rsamtools::indexBam(bam_path)
        print("baifile already generated.")
    }

    ## build dataframe of read bin
    # annotate bin
    which <- getAnnotate_bin(bin_size, ref_genome, blacklist)
    if(check_Chr_UCSU_format(bam_path)){
        which <- convert_GRch2UCSC(which)
    }

    # flag
    flag <- Rsamtools::scanBamFlag(
        isPaired = TRUE,
        isUnmappedQuery = FALSE,
        isDuplicate = FALSE,
        isMinusStrand = FALSE,
        hasUnmappedMate = FALSE,
        isSecondaryAlignment = FALSE,
        isMateMinusStrand = TRUE
    )

    # parameters
    param <- Rsamtools::ScanBamParam(
        what = c("rname","pos","cigar","isize","seq"),
        flag = flag,
        which = which,
        mapqFilter = min_mapQ
    )

    # merge all parts
    print("Reading bam file.")
    bam <- Rsamtools::scanBam(
        file = bam_path,
        index = bam_path,
        param = param
    )

    # extract read on target bin
    if(!is.null(target_bed)) {
        print("Extracting read on target segments.")
        bam <- extract_read_on_target(bam, target_bed)
    }

    class(bam)="BAM"

    if(check_Chr_UCSU_format(bam_path)){
        bam <- convert_UCSC2Chr(bam)
    }

    return(bam)
}

extract_read_on_target <- function(bam, target_bed) {
    if(check_gzip(target_bed)){
        target_df <- as.data.frame(read.table(gzfile(target_bed),
        header = FALSE, 
        stringsAsFactors = FALSE,
        sep = "\t"))[,seq_len(3)]
    } else {
        target_df <- as.data.frame(read.table(target_bed,
        header = FALSE, 
        stringsAsFactors = FALSE,
        sep = "\t"))[,seq_len(3)]
    }

    target_gr=GenomicRanges::GRanges(
    seqnames = target_df[, 1],
    ranges = IRanges::IRanges(start = as.numeric(target_df[, 2]), end=as.numeric(target_df[, 3]))
    )

    read_on_target_bam = lapply(bam, function(region){
        bin_gr = GenomicRanges::GRanges(
            seqnames = as.character(region$rname),
            ranges = IRanges::IRanges(start = region$pos, end = region$pos + region$qwidth),
            qname = region$qname,
            rname = region$rname,
            pos = region$pos,
            isize = region$isize,
            seq = region$seq)
        
        read_on_target_gr <- GenomicRanges::findOverlaps(bin_gr, target_gr)

        if(length(read_on_target_gr@from) == 0 )
            read_on_target_bam_gr = bin_gr
        else{
            read_on_target_bam_gr = bin_gr[read_on_target_gr@from]
        }

        return_lst = list("qname" = read_on_target_bam_gr$qname,
            "rname" = read_on_target_bam_gr$rname,
            "pos" = read_on_target_bam_gr$pos,
            "isize" = read_on_target_bam_gr$isize,
            "seq" = read_on_target_bam_gr$seq)
    })
    return(read_on_target_bam)
}
