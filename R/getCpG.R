#' @export

getCpG <- function(read_bin) {
    total <- length(read_bin)
    range <- 1:total
    count <- list()
    cat("Total of segment bin: ", total, "\n")

    for(i in range) {
        start <- Sys.time()
        d <- read_bin[[i]][[4]]
        nn <- triNUC_count(d)
        count[[i]] <- nn
        cat("Step ", i, ":\n")
        print(count[[i]])
        print( Sys.time() - start )
    }

    new_CG <- Reform(count)

    return(new_CG)
}

triNUC_count <- function(bin) {

    cg <- length(grep("^CG.*", bin, value=TRUE))
    ag <- length(grep("^AG.*", bin, value=TRUE))
    tg <- length(grep("^TG.*", bin, value=TRUE))
    gg <- length(grep("^GG.*", bin, value=TRUE))
    cc <- length(grep("^CC.*", bin, value=TRUE))
    ac <- length(grep("^AC.*", bin, value=TRUE))
    tc <- length(grep("^TC.*", bin, value=TRUE))
    gc <- length(grep("^GC.*", bin, value=TRUE))
    ca <- length(grep("^CA.*", bin, value=TRUE))
    aa <- length(grep("^AA.*", bin, value=TRUE))
    ta <- length(grep("^TA.*", bin, value=TRUE))
    ga <- length(grep("^GA.*", bin, value=TRUE))
    ct <- length(grep("^CT.*", bin, value=TRUE))
    at <- length(grep("^AT.*", bin, value=TRUE))
    tt <- length(grep("^TT.*", bin, value=TRUE))
    gt <- length(grep("^GT.*", bin, value=TRUE))

    acg <- length(grep("^ACG.*", bin, value=TRUE))
    tcg <- length(grep("^TCG.*", bin, value=TRUE))
    ccg <- length(grep("^CCG.*", bin, value=TRUE))
    gcg <- length(grep("^GCG.*", bin, value=TRUE))


    all <- list(cg, ag, tg, gg, cc, ac, tc, gc, ca, aa, ta, ga, ct, at, tt, gt, acg, tcg, ccg, gcg)

    # CpG_ratio <- cgn/ncg

    return(all)
}

CpG_ratio <- function(list_tri) {
    cg <- list_tri[[1]][1]
    acg <- list_tri[[17]][1]
    tcg <- list_tri[[18]][1]
    ccg <- list_tri[[19]][1]
    gcg <- list_tri[[20]][1]

    CpG <- cg/(acg + tcg + ccg + gcg)

    return(CpG)
}

Reform <- function(CpG_list) {
    mat <- matrix(unlist(CpG_list), ncol = 1, byrow = TRUE)
    df <- as.data.frame(mat)
    names(df) <- c("CpG_ratio")

    return(df)
}