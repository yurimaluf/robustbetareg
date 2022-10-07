# .onLoad <- function(...){
#   #msg <- paste("Loading", pkgname)
#   cat("Em caso de algum bug por favor reportar \n")
# }

.onAttach<-function(...){
  packageStartupMessage("In case of a bug please report it to: yurimaluf@gmail.com \n")
  #cat("In case of a bug please report it to: yurimaluf@gmail.com \n")
  #rstudioapi::showDialog( "BugReports","Report bugs at:", url="https://github.com/yurimaluf/RobustBetaReg/issues")
}
