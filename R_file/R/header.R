.onAttach <- function(...) {

   mydate <- date()
   x <- regexpr("[0-9]{4}", mydate)
   this.year <- substr(mydate, x[1], x[1] + attr(x, "match.length") - 1)

   packageStartupMessage("\n## Ordered Optimal Classification: Ideal Point Estimation Package")
   packageStartupMessage("## Christopher Hare and Keith T. Poole")

}

.onUnload <- function(libpath) {
    library.dynam.unload("ooc", libpath)
}
