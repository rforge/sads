## Correcoes na plot.count do untb para deixar os rotulos das especies ok

## como colocar lables obliquos em x, faq 7.27
 ## Increase bottom margin to make room for rotated labels
     par(mar = c(7, 4, 4, 2) + 0.1)
     ## Create plot with no x axis and no x axis label
     plot(1 : 8, xaxt = "n",  xlab = "")
     ## Set up x axis with tick marks alone
     axis(1, labels = FALSE)
     ## Create some text labels
     labels <- paste("Label", 1:8, sep = " ")
     ## Plot x axis labels at default tick marks
     text(1:8, par("usr")[3] - 0.25, srt = 45, adj = 1,
          labels = labels, xpd = TRUE)
     ## Plot x axis label at line 6 (of 7)
     mtext(1, text = "X Axis Label", line = 6)

##tentando implementar no plot.count

require(untb)

 
plot.count <- function (x, uncertainty = FALSE, expectation = FALSE, theta = NULL, n = 10, ...) {
  x <- as.count(x)
  ## Increase bottom margin to make room for rotated labels
  par(mar = c(7, 4, 4, 2) + 0.1)
  ## Create plot with no x axis and no x axis label
  plot.default(x, log = "y", col = "red", pch = 16, type = "b", 
               xaxt="n",ylab = "abundance", xlab="", ...)
  ## Set up x axis with tick marks alone
  ind <- seq.int(1,length(x),length.out=min(length(x),10))
  axis(1, labels = FALSE, at=ind)
  axis(2)
  ## Plot x axis labels at default tick marks
  text(ind, 10^(par("usr")[3] - 0.25), srt = 30, adj = 1,
       labels = names(x)[ind], xpd = TRUE)
  ## Plot x axis label at line 6 (of 7)
  mtext(1, text = "Species Rank", line = 6)
  if (uncertainty | expectation) {
    J <- no.of.ind(x)
    if (is.null(theta)) {
      theta <- optimal.theta(x)
    }
  }
  if (uncertainty) {
    for (i in 1:n) {
      jj <- rand.neutral(J = J, theta = theta)
      points(1:no.of.spp(jj), jj, type = "l", col = "gray")
    }
  }
  if (expectation) {
    points(1:J, expected.abundance(J = J, theta = theta), 
           type = "b", pch = 16)
  }
}

data(saunders)
.GlobalEnv:plot.count(extant(extractor(saunders,1)))
