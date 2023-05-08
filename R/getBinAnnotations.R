#' @export 

getAnnotate_bin <- function(bin_size, ref_genome, blacklist) {
    # seperte bin in each size
    if(ref_genome == 'hg38') {
        bins <- QDNAseq::getBinAnnotations(binSize=bin_size, genome=ref_genome)
    } else {
        stop("Error: Your type of reference genome is valid.")
    }

    # build bin dataframe and filter blacklist out
    get_bins <- as.data.frame(bins@data)
    if(blacklist == TRUE) {
        get_bins <- get_bins[which(get_bins$chromosome!="Y" & get_bins$mappability>=1 & get_bins$blacklist==0),]
    } else if(blacklist == FALSE) {
        get_bins <- get_bins[which(get_bins$chromosome!="Y" & get_bins$mappability>=1),]
    }
    get_bins_gr <- GRanges(
        seqnames = get_bins$chromosome,
        ranges = IRanges(start= get_bins$start, end = get_bins$end),
        gc = get_bins$gc,
        mappability = get_bins$mappability
    )
    return(get_bins_gr)
}