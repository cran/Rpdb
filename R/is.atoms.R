is.atoms <- function(x)
{
  to.return <- any(attr(x,which="class") == "atoms")
  return(to.return)
}