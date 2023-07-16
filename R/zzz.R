.onAttach <- function(libname, pkgname) {
  title <- read.dcf(file.path(libname, pkgname, "DESCRIPTION"), "Title")
  title <- as.character(title)
  packageStartupMessage(title, " : Version ", utils::packageVersion("fPASS"))
}

