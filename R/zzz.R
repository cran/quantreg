# Added function:

".First.lib" <- 
function(lib, pkg) {
   library.dynam("quantreg", pkg, lib) }
