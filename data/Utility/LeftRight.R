library(GenomicAlignments)

# in newer package versions, left() and right() are deprecated/defunct and are
# replaced by first() and last
if ( packageVersion("GenomicAlignments") >= "1.6.0") {

# these are the orignial functions from GenomicAlignments, taken from:
# https://github.com/Bioconductor-mirror/GenomicAlignments/blob/release-3.0/R/utils.R

    invertRleStrand <- function(x)
    {
        x_strand <- strand(x)
        runValue(x_strand) <- strand(runValue(x_strand) == "+")
        strand(x) <- x_strand
        x
    }



    right <- function(x){
        x_first <- x@first
        x_last <- invertRleStrand(x@last)
    
        right_is_first <- which(strand(x_first) == "-")
        idx <- seq_len(length(x))
        idx[right_is_first] <- idx[right_is_first] + length(x)
    
        ans <- c(x_last, x_first)[idx]
        setNames(ans, names(x))
    }


    left <- function(x){
        x_first <- x@first
        x_last <- invertRleStrand(x@last)
        
        left_is_last <- which(strand(x_first) == "-")
        idx <- seq_len(length(x))
        idx[left_is_last] <- idx[left_is_last] + length(x)
    
        ans <- c(x_first, x_last)[idx]
        setNames(ans, names(x))
    }


}
