f2g <- function(mypdf, mypng) { # convert pdf2png
  mypng <- sub('.png', '', mypng)
  system(paste0('/usr/bin/pdftoppm -png ', mypdf, ' ', mypng, ' -singlefile'))
}
