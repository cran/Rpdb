is.coords <- function(x)
{
  to.return <- any(class(x) == "coords")
  return(to.return)
}