#' A function to read a HOMER-formatted peak file into a GRanges object
#' 
#' A function to read a HOMER formatted peak/position file into a GRanges object. HOMER is a popular DNA motif analysis and general NGS program suite.
#' Peak/position files have a commented header and are formatted as:
#' PeakID, chr, start, end, strand ...
#' @param file an absolute or relative path to a HOMER-formatted peak file
#' @param header.regex an optional regular expression to identify the header line (Default:"^#PeakID").
#' @param meta.cols \code{named list} that maps column numbers to meta data columns. e.g. list(name=5, score=10), which means 5th column will be named "name", and 10th column will be named "score" and their contents will be a part of the returned \code{GRanges object}.
#' @param keep.all.metadata \code{logical} determining if the extra columns ( the ones that are not designated by chr,start,end,strand and meta.cols arguments ) should be kept or not. (Default:\code{FALSE})
#' @param remove.unusual if \code{TRUE}(default) remove the chromosomes with unsual names, such as chrX_random. (Default:\code{FALSE})
#' @return \code{\link{Granges}} object.
#' @seealso \code{\link{readGeneric}}
#' @examples
#' # peaks <- readHomer("peaks.txt", keep.all.metadata=TRUE)
#' # head(peaks)
readHomer <- function(file, header.regex="^#PeakID", meta.cols = NULL, keep.all.metadata = FALSE, remove.unusual = FALSE)
{
  if (keep.all.metadata) {
    # Because header is commented, we need to regex it and pass
    # the identifiers to mcols
    track <- scan(file, n=50, what="character", sep="\n", quiet=T)
    homer.header <- track[grepl(header.regex, track)]
    homer.header <- sub("#", "", homer.header)
    mcols <- strsplit(homer.header, split = "\t")
    colNums <- as.vector(1:length(unlist(mcols)))
    names(colNums) <- unlist(mcols)
    mcols <- colNums
    mcols <- mcols[-2] # Delete chr as not caught by blacklist
  } else { mcols <- meta.cols}
  
  g <- genomation::readGeneric(file, chr=2, start = 3, end = 4, strand = 5, meta.cols = mcols, remove.unusual = remove.unusual)
  return(g)
}