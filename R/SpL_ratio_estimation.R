#' @export

SpL_calculation <- function(
    start <- Sys.time()
    read_bin,
    ref_genome = "hg38",
    short = c(90,150),
    long = c(151,250)
) {
    total <- length(read_bin)
    range <- 1:total
    SpL_ratio <- list()
    pos <- Position_extract(read_bin)

    for(i in range) {
        bin <- read_bin[[i]][[3]]
        chr_bin <- droplevels(read_bin[[i]][[1]])
        chr_1st_inBin <- toString(read_bin[[i]][[1]][1])
        chr_ext <- CHR_read(chr_bin, chr_1st_inBin)
        short_range <- length(which(bin >= short[1] & bin <= short[2]))
        long_range <- length(which(bin > long[1] & bin <= long[2]))
        SpL <- SpL_func(short_range, long_range)
        SpL_ratio[[i]] <- list(chr_ext, pos[i], SpL)
    }

    new_form <- Reform_data(SpL_ratio)
    print( Sys.time() - start )
    
    return(new_form)
}

SpL_func <- function(short, long) {
	ratio = short / long
	if(is.nan(ratio)) {
	return('-')
	} else {
	return(ratio)
	}
}

CHR_read <- function(input, chr_no) {
    check_chr <- all(input %in% list(chr_no))
    if(check_chr == TRUE) {
        return(chr_no)
    } else if(check_chr == FALSE) {
        return('-')
    } else if(is.na(check_chr)) {
        return('-')
    }
}

Reform_data <- function(spl_list) {
    mat <- matrix(unlist(spl_list), ncol = 3, byrow = TRUE)
    df <- as.data.frame(mat)
    names(df) <- c("chromosome", "position", "short_long")

    return(df)
}

Position_extract <- function(read) {
    pos <- names(read)
    pos_del <- gsub(".*:","",pos)
    return(pos_del)
}