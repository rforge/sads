octav.rel <- function (x, oct, ...) 
{
    if (missing(oct)) {
        oct <- 1:(ceiling(max(log2(x))) + 1)
        if (any(x < 1)) {
            octlower <- ceiling(min(log2((x))) + 1):0
            oct <- c(octlower, oct)
        }
    }
    n <- 2^(oct - 1)
    oc.class <- cut(x, breaks = c(0, n), labels = oct)
    res <- as.data.frame(table(oc.class))
    res$upper <- n
    names(res)[1] <- "octave"
    res$Freq.rel <- res[,2]/sum(res[,2])
    res[, c(1, 3, 2, 4)]
}
