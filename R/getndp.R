`getndp` <-
function(x, tol=2*.Machine$double.eps)
{
  ndp <- 0
  while(!isTRUE(all.equal(x, round(x, ndp), tol=tol))) ndp <- ndp+1
  if(ndp > -log10(tol)) warning("Tolerance reached, ndp possibly
underestimated.")
  ndp
}

