library(DensityBAST)
library(spatstat)
library(doRNG)

data(chicago)
plot(chicago)

chicago_intensity_function <- function(x, y,L=chicago$domain) {
  # Create a temporary point
  pt <- ppp(x, y, window = Frame(L))  # spatial point (x,y)
  chicago_intens=apply(getLinearDensitySamples(chicago.fit,pt),2,mean)
}

chicago_intensity_function_sd <- function(x, y,L=chicago$domain) {
  # Create a temporary point
  pt <- ppp(x, y, window = Frame(L))  # spatial point (x,y)
  chicago_intens=apply(getLinearDensitySamples(chicago.fit,pt),2,sd)
}

chicago_linim <- linim(chicago$domain, chicago_intens, locations = tesselation)
## fit chicago data using densitybast
chicago.fit=DensityBASTFit_lpp(chicago)
#
Z <- as.im(chicago_intensity_function, Frame(chicago$domain))
intensity.mean <- linim(chicago$domain, Z)
plot(intensity.mean)

Z <- as.im(chicago_intensity_function_sd, Frame(chicago$domain))
intensity.sd <- linim(chicago$domain, Z)
plot(intensity.sd)
