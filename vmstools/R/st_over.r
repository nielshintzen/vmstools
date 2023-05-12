st_over <- function(x, y) {
  require(sf)
  sapply(sf::st_intersects(x, y), function(z)
    if (length(z) == 0) {
      NA_integer_
    } else {
      z[1]
    })
}