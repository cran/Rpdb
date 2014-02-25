is.cryst1 <- function(x)
{
  to.return <- any(attr(x,which="class") == "cryst1")
  return(to.return)
}