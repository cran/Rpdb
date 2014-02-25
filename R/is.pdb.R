is.pdb <- function(x)
{
  to.return <- any(class(x) == "pdb")
  return(to.return)
}