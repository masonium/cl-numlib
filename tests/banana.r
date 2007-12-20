banana <- function(x) {
  -100*(x[2]-x[1]^2)^2-(1-x[1])^2
}

bananagr <- function(x,y) {
  l <- x[2]-x[1]^2
  c(2*(1-x[1])+400*l*x[1], -200*l)
}

optim(c(1,1), banana, bananagr, method="BFGS",
      control=list(fnscale=-1,trace=3))
