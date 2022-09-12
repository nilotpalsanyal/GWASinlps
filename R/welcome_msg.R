.onAttach <- function(libname, pkgname) {
  message <- c("\n Welcome to GWASinlps...Happy selection!",
               "\n \n Website: https://nilotpalsanyal.github.io/GWASinlps/",
               "\n Bug report: https://github.com/nilotpalsanyal/GWASinlps/issues")
  packageStartupMessage(message)
}