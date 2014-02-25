is.conect <- function(x)
{
  to.return <- any(attr(x,which="class") == "conect")
  return(to.return)
}