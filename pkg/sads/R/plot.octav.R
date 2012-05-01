plot.octav <- function(x,...){
  barplot(height=x$Freq, space=0, ylab="Frequency", xlab="Abundance class",...)
  axis(1, at=x$octave, labels=x$upper)
}
