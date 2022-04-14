.onLoad <- function (libname, packagename)
{
  # these are all the data files in your package
  files <- c(
    "trueLabel.RData",
    "sig.RData",
    "labels.RData",
    "dataSC_1.Rdata",
    "dataSC_2.RData",
    "dataBulk.RData"
  )
  package.folder <- file.path("extdata")
  network.folder <- file.path("inst", "extdata")


  package.files <- file.path(package.folder, files)
  network.files <- file.path(network.folder, files)


  # now we need to check that 'package.files' actually exist
  # start by finding the package on the computer
  #
  # in folder 'libname'
  # look for package 'pkgname'
  package.location <- system.file(lib.loc = libname, package =
                                    pkgname, mustWork = TRUE)


  # the location of DWLS on github (or maybe some other file storage website)
  # if on github, make sure to use raw.githubusercontent.com and NOT
  #github.com
  network.location <-
    "https://github.com/sistia01/DWLS/main"


  package.paths <- file.path(package.location, package.files)
  network.paths <- file.path(network.location, network.files)


  # which of 'package.paths' DO NOT EXIST?
  # these files which don't exist are the files that we need to download
  need2download <- !file.exists(package.paths)
  if (any(need2download)) {


    # get the files which we need to download
    package.paths <- package.paths[need2download]
    network.paths <- network.paths[need2download]


    # one more unfortunate thing. on Windows, the download mode
    #will auto-correct
    # to `mode = "w"` or `mode = "wb"` when necessary, but not on
    #Unix-alikes
    # here is the solution, copied directly from R source
    is.binary <- grepl("\\.(gz|bz2|xz|tgz|zip|jar|rd[as]|RData)$",
                       network.paths)
    modes <- ifelse(is.binary, "wb", "w")


    # just one more thing before you actually download the files.
    # we need to make sure that the directories where we will be
    # placing these
    # files actually exists:
    folders <- unique(dirname(package.paths))
    for (folder in folders)
      dir.create(folder, showWarnings = FALSE, recursive = TRUE)


    # if you haven't seen .mapply before,
    # it's useful for looping over multiple variables at once
    # you give it a data.frame, and it will loop over each row,
    # passing as arguments to a function
    arguments <- data.frame(
      url      = network.paths,
      destfile = package.paths,
      quiet    = TRUE,
      mode     = modes
    )
    .mapply(utils::download.file, arguments, NULL)
    # .mapply will do something along the lines of:
    #
    # utils::download.file(
    #     url      = arguments[1, "url"],
    #     destfile = arguments[1, "destfile"],
    #     quiet    = arguments[1, "quiet"],
    #     mode     = arguments[1, "mode"]
    # )
    # utils::download.file(
    #     url      = arguments[2, "url"],
    #     destfile = arguments[2, "destfile"],
    #     quiet    = arguments[2, "quiet"],
    #     mode     = arguments[2, "mode"]
    # )
    #
    # until it reaches the end of the data.frame
  }
}
