#' A function to load some data when the package opens.
#'
#' @param libname Library name.
#' @param pkgname Package name.
#' @export

.onAttach = function(libname=NULL, pkgname="STRAND") {
  github_auth_token = function() {
  github_pat = Sys.getenv("GITHUB_PAT")
  if (nzchar(github_pat)) {
    auth_token = c(Authorization = paste0("token ", github_pat))
  } else {
    auth_token = NULL
  }
  auth_token
}

try_download = function(download_url, destination_file,
                          quiet = TRUE) {
  download_status = try(
    suppressWarnings(
      utils::download.file(url = download_url,
                           destfile = destination_file,
                           quiet = quiet,
                           headers = github_auth_token())
    ),
    silent = TRUE
  )
  download_status
}

download_with_retries = function(download_url,
                                  destination_file,
                                  retries = 5,
                                  pause_sec = 5,
                                  quiet = TRUE) {
    download_rc = try_download(download_url, destination_file,
                                quiet = quiet)
    num_retries = 0
    while (num_retries < retries && inherits(download_rc, "try-error")) {
      Sys.sleep(pause_sec)
      num_retries = num_retries + 1
      download_rc = try_download(download_url, destination_file, quiet = quiet)
    }
    download_rc
}

latest_released_version = function(quiet=TRUE, ...) {
  dest_file = tempfile(pattern = "releases-", fileext = ".json")
  download_url = "https://api.github.com/repos/ctross/STRAND/releases/latest"
  release_list_downloaded = download_with_retries(download_url, dest_file, quiet = quiet, ...)
  if (inherits(release_list_downloaded, "try-error")) {
    stop("GitHub download of release list failed with error: ",
        attr(release_list_downloaded, "condition")$message,
        call. = FALSE)
  }
  release = jsonlite::read_json(dest_file)
  sub("v", "", release$tag_name)
}

  packageStartupMessage("Bei diesem Spaziergang an den STRAND scharfen wir unsere Sinne fur die Sternbilder hoch am Himmel!")

  if(pingr::is_online()){
  latest_version = latest_released_version(retries = 0)

  installed_version = "the_height_of_our_halo"
  
  packageStartupMessage(paste0("This version of STRAND is: ", installed_version, "."))

  if(installed_version != latest_version){
  packageStartupMessage(paste0("The latest release of STRAND is: ", latest_version, ".\n 
   STRAND is under active development, so consider updating with: devtools::install_github('ctross/STRAND@", latest_version,"')"))
  }}
}







