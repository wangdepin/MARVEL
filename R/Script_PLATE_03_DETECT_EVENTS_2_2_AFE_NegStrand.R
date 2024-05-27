#' @title Detect alternative first exons on negative strand
#'
#' @description Detects alternative first exons, specifically for genes transcribed on the negative strand of the DNA.
#'
#' @param MarvelObject S3 object generated from \code{CreateMarvelObject} function.
#' @param parsed.gtf Data frame. GTF file with the gene_id parsed. Generated from the \code{DetectEvents.AFE} function.
#' @param min.cells Numeric value. The minimum number of cells in which the gene is expressed for the gene to be included for splicing event detected and quantification. To be used in conjunction with \code{min.expr} argument. Default value is \code{50}.
#' @param min.expr Numeric value. The minimum expression value for the gene to be considered to be expressed in a cell. Default value is \code{1}.
#' @param track.progress Logical. If set to \code{TRUE}, progress bar will appear to track the progress of the rate-limiting step of this function, which is the extraction of the final exon-exon junctions. Default value is \code{FALSE}.
#'
#' @return An object of class S3 with new slot \code{MarvelObject$SpliceFeature$AFE.NegStrand}.
#'
#' @importFrom plyr join
#' @import methods
#' @importFrom utils txtProgressBar setTxtProgressBar
#'
#' @export
#'
#' @examples
#' marvel.demo <- readRDS(system.file("extdata/data", "marvel.demo.rds", package="MARVEL"))
#'
#' marvel.demo <- DetectEvents.AFE.NegStrand(MarvelObject=marvel.demo,
#'                                           parsed.gtf=NULL,
#'                                           min.cells=5,
#'                                           min.expr=1,
#'                                           track.progress=FALSE
#'                                           )

DetectEvents.AFE.NegStrand <- function(MarvelObject, parsed.gtf=NULL, min.cells=50, min.expr=1, track.progress=FALSE) {

    # Define arguments
    df <- parsed.gtf
    df.sj <- MarvelObject$SpliceJunction
    df.gene <- MarvelObject$Exp
    min.cells <- min.cells
    min.expr <- min.expr
    track.progress <- track.progress
    
    # Check if GTF provided
    if(is.null(df)) {
        message("Please provide GTF into the parsed.gtf option")
        return(MarvelObject)
    }
    
    # Create row names
    row.names(df.sj) <- df.sj$coord.intron
    df.sj$coord.intron <- NULL
    row.names(df.gene) <- df.gene$gene_id
    df.gene$gene_id <- NULL
    
    # Recode SJ NA's as 0'
    df.sj[is.na(df.sj)] <- 0
    
    # Retrieve gene_id metadata for annotation later
    message("Retrieving gene metadata from GTF...")
    df.small <- df[which(df$V3=="gene"), ]
    attr <- strsplit(df.small$V9, split=";")
    df.small$gene_id <- textclean::mgsub(sapply(attr, function(x) grep("gene_id", x, value=TRUE)), c("gene_id", " ", "\""), "")
    df.small$gene_short_name <- textclean::mgsub(sapply(attr, function(x) grep("gene_name", x, value=TRUE)), c("gene_name", " ", "\""), "")
    df.small$gene_type <- textclean::mgsub(sapply(attr, function(x) grep("gene_biotype", x, value=TRUE)), c("gene_biotype", " ", "\""), "")
    df.feature <- df.small[, c("gene_id", "gene_short_name", "gene_type")]
    
    ##########################################################################
    ######################### SUBSET EXPRESSED GENES #########################
    ##########################################################################

    message("Retrieving expressed genes...")

    gene_ids <- names(apply(df.gene, 1, function(x) {sum(x >= min.expr)}))[which(apply(df.gene, 1, function(x) {sum(x >= min.expr)}) >= min.cells)]
    df <- df[which(df$V7=="-"), ]
    gene_ids.overlap <- intersect(gene_ids, unique(df$gene_id))
    df <- df[which(df$gene_id %in% gene_ids.overlap), ]

    message(paste(length(gene_ids.overlap), " expressed genes identified", sep=""))
    
    ##########################################################################
    ########################### SUBSET FIRST SJ ##############################
    ##########################################################################

    df <- df[which(df$V3=="exon"), ]
    attr <- strsplit(df$V9, split=";")
    df$transcript_id <- textclean::mgsub(sapply(attr, function(x) grep("transcript_id", x, value=TRUE)), c("transcript_id", " ", "\""), "")
    
    freq <- as.data.frame(table(df$transcript_id))
    transcript_ids <- freq[which(freq$Freq >=4), 1]
    df <- df[which(df$transcript_id %in% transcript_ids), ]
    
    message(paste("Retrieving final exon-exon junctions from ", length(transcript_ids), " multi-exon transcripts", sep=""))
    transcript_ids <- unique(df$transcript_id)
    .list <- list()
    if(track.progress==TRUE) {
        pb <- txtProgressBar(1, length(transcript_ids), style=3)
    }
    for(i in 1:length(transcript_ids)) {
        df.small <- df[which(df$transcript_id %in% transcript_ids[i]), ]
        df.small <- df.small[order(df.small$V4, decreasing=TRUE), ]
        if (nrow(df.small) > 1) {
            .list[[i]] <- data.frame("gene_id"=df.small$gene_id[1],
                                     "chr"=df.small$V1[1],
                                     "V1"=df.small$V4[2],
                                     "V2"=df.small$V5[2],
                                     "V3"=df.small$V4[1],
                                     "V4"=df.small$V5[1],
                                     stringsAsFactors=FALSE
            )
        }
        if(track.progress==TRUE) {
            setTxtProgressBar(pb, i)
        }
    }
    df <- do.call(rbind.data.frame, .list)
    df <- unique(df)
    
    par_y <- grep("_PAR_Y", df$gene_id)
    if(length(par_y) != 0) {
        df <- df[-grep("_PAR_Y", df$gene_id), ]
    }
    
    df$coord.intron <- paste(df$chr, df$V2 + 1, df$V3 - 1, sep=":")
    df <- df[which(df$coord.intron %in% row.names(df.sj)), ]

    ##########################################################################
    ############################# REMOVE A5SS ################################
    ##########################################################################

    coords <- unique(df$V4)
    coords.no.a5ss <- NULL
    for(i in 1:length(coords)) {
        df.small <- df[which(df$V4==coords[i]), ]
        if(length(unique(df.small$V3))==1) {
            coords.no.a5ss[i] <- coords[i]
        }
    }
    df <- df[which(df$V4 %in% coords.no.a5ss), ]
    
    ##########################################################################
    ######################### COLLAPSE COORDINATES ###########################
    ##########################################################################

    message("Collapsing redundant coordinates/exons...")
    coord.introns <- unique(df$coord.intron)
    .list <- list()
    for(i in 1:length(coord.introns)) {
        df.small <- df[which(df$coord.intron %in% coord.introns[i]), ]
        if(nrow(df.small) == 0) next
        if(nrow(df.small)==1) {
            .list[[i]] <- df.small
        } else {
            if (!all(is.na(df.small$V1))) {
                df.small$V1 <- max(df.small$V1, na.rm=TRUE)
            }
            if (!all(is.na(df.small$V4))) {
                df.small$V4 <- min(df.small$V4, na.rm=TRUE)
            }
            .list[[i]] <- df.small[1, ]
        }
    }
    df <- do.call(rbind.data.frame, .list)
    
    freq <- as.data.frame(table(df$gene_id))
    names(freq) <- c("gene_id", "n.transcripts")
    df <- join(df, freq, by="gene_id", type="left")
    df <- df[which(df$n.transcripts >= 2), ]
    df$n.transcripts <- NULL
    
    coords <- unique(df$V2)
    .list <- list()
    for(j in 1:length(coords)) {
        df.small <- df[which(df$V2 %in% coords[j]), ]
        if(nrow(df.small)==1) {
            .list[[j]] <- NULL
        } else if(nrow(df.small) == 2) {
            df.small <- df.small[order(df.small$V3, decreasing=TRUE), ]
            tran_id <- paste(df.small$chr[1], ":", df.small$V3[1], ":", df.small$V4[1],
                             "|",
                             df.small$V3[2], ":", df.small$V4[2],
                             ":-@",
                             df.small$chr[1], ":", df.small$V1[1], ":", df.small$V2[1],
                             sep=""
            )
            .list[[j]] <- data.frame("tran_id"=tran_id,
                                     "gene_id"=df.small$gene_id[1],
                                     stringsAsFactors=FALSE
            )
        } else if(nrow(df.small) >=3) {
            df.small <- df.small[order(df.small$V3, decreasing=TRUE), ]
            df.small. <- df.small[-1, ]
            tran_ids <- NULL
            for(i in 1:nrow(df.small.)) {
                df.small.. <- df.small.[i, ]
                tran_ids[i] <- paste(df.small$chr[1], ":", df.small$V3[1], ":", df.small$V4[1],
                                     "|",
                                     df.small..$V3[1], ":", df.small..$V4[1],
                                     ":-@",
                                     df.small$chr[1], ":", df.small$V1[1], ":", df.small$V2[1],
                                     sep=""
                )
            }
            .list[[j]] <- data.frame("tran_id"=tran_ids,
                                     "gene_id"=df.small$gene_id[1],
                                     stringsAsFactors=FALSE
            )
        }
    }
    df <- do.call(rbind.data.frame, .list)
    
    df.feature <- join(df, df.feature, by="gene_id", type="left")
    message(paste(nrow(df.feature), " AFE identified", sep=""))
    
    df.feature$event_type <- "AFE"
    col.others <- names(df.feature)[-which(names(df.feature) %in% c("tran_id", "event_type"))]
    df.feature <- df.feature[, c("tran_id", "event_type", col.others)]
    df.feature <- unique(df.feature)
    MarvelObject$SpliceFeature$AFE.NegStrand <- df.feature
    return(MarvelObject)
}
