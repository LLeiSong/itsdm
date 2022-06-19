.onAttach <- function(libname, pkgname) {
  pkg.version <- utils::packageVersion("itsdm")
  # msg <- c(
  #   "  _    _         \n",
  #   " / \\__/ \\_____     _ _           _\n",
  #   "/  /  \\  \\    `\\  (_) |_ ___  __| |_ __ ___\n",
  #   ")  \\''/  (     |\\ | | __/ __|/ _` | '_ ` _ \\\n",
  #   "`\\__)/__/'_\\  / ` | | |_\\__ \\ (_| | | | | | |\n",
  #   sprintf("   //_|_|~|_|_|   |_|\\__|___/\\__,_|_| |_| |_| Version %s\n",
  #           pkg.version))

  msg <- c(
    " _ _           _\n",
    "(_) |_ ___  __| |_ __ ___\n",
    "| | __/ __|/ _` | '_ ` _ \\\n",
    "| | |_\\__ \\ (_| | | | | | |\n",
    sprintf("|_|\\__|___/\\__,_|_| |_| |_| Version %s\n",
            pkg.version))
  packageStartupMessage(msg)
}
