##load data
#either normalized or not
path <- "/home/s3/afr/data/human/E2F7/genome_coverage/"
outpath <- "/home/s3/afr/data/human/E2F7/"
dataTable.norm <- read.delim(paste(path, "genome_coverage_norm.bed", sep=""), sep="\t", header=TRUE)


panel.hist <- function(x, ...) {
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col="gray", ...)
}
panel.cor <- function(x, y, digits=2, prefix="", cex.cor) {
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- cor(x, y, use="pairwise.complete")
  txt <- format(c(r, 0.123456789), digits=digits)[1]
  txt2 <- format(c(abs(r), 0.123456789), digits=digits)[1]
  txt <- paste(prefix, txt, sep="")
  if(missing(cex.cor)) cex <- 0.8/strwidth(txt2)
  text(0.5, 0.5, txt, cex = cex * abs(r))
}
mypoints <- function(x, y, n=3000, ...) {
  ok <- which(is.finite(x) & is.finite(y))
  if(length(ok)>n) ok <- sample(ok, n)
  x <- x[ok]; y <- y[ok]
  points(x, y, pch=".")
}

mypoints2 <- function(...) {par(new=TRUE);smoothScatter(..., nrpoints=0)}

png("ChIP.png", height=2000, width=2500, pointsize=20)
pairs(dataTable.norm[,c(-1,-2,-3)], diag.panel=panel.hist, upper.panel=panel.cor, lower.panel=mypoints2)
dev.off()


# correlations at gene level