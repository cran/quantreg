# Added function:

".First.lib" <-
function(lib, pkg) {
   library.dynam("quantreg", pkg, lib)
   print("quantreg library 3.01loaded")}
