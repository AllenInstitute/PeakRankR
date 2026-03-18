#' Check That bedtools Is Accessible
#'
#' Verifies that the `bedtools` binary can be found, either on the system
#' `PATH` or at the path set via `options(bedtools.path = ...)`. If bedtools
#' cannot be found, the function stops with a clear, actionable error message.
#'
#' @param quiet Logical. If `TRUE`, suppresses the success message.
#'   Default: `FALSE`.
#'
#' @return Invisibly returns the resolved path to bedtools if found.
#'
#' @details
#' You can point PeakRankR to a specific bedtools installation by setting:
#' ```r
#' options(bedtools.path = "/path/to/bedtools/bin/")
#' ```
#' If this option is not set, the function uses `Sys.which("bedtools")` to
#' locate the binary on your `PATH`.
#'
#' To install bedtools, see:
#' <https://bedtools.readthedocs.io/en/latest/content/installation.html>
#'
#' @examples
#' \dontrun{
#' # Auto-detect bedtools on PATH
#' check_bedtools()
#'
#' # Specify a custom path
#' options(bedtools.path = "/usr/local/bin/")
#' check_bedtools()
#' }
#'
#' @export
check_bedtools <- function(quiet = FALSE) {

  # Check for user-specified path
  opt_path <- getOption("bedtools.path", default = NULL)

  if (!is.null(opt_path)) {
    candidate <- file.path(opt_path, "bedtools")
    # On Windows, also try .exe
    if (.Platform$OS.type == "windows") {
      candidate_win <- paste0(candidate, ".exe")
      if (file.exists(candidate_win)) {
        if (!quiet) message("bedtools found at: ", candidate_win)
        return(invisible(candidate_win))
      }
    }
    if (file.exists(candidate)) {
      if (!quiet) message("bedtools found at: ", candidate)
      return(invisible(candidate))
    }
    warning(
      "bedtools not found at the path set in options(bedtools.path): ",
      opt_path,
      "\nFalling back to PATH search."
    )
  }

  # Try system PATH
  bt_path <- Sys.which("bedtools")

  if (nchar(bt_path) == 0 || bt_path == "") {
    stop(
      "bedtools was not found on your system PATH.\n\n",
      "To fix this, either:\n",
      "  1. Install bedtools: https://bedtools.readthedocs.io/en/latest/content/installation.html\n",
      "  2. Set its location with: options(bedtools.path = \"/path/to/bedtools/bin/\")\n\n",
      "To find an existing installation, run in your terminal:\n",
      "  which bedtools     # macOS / Linux\n",
      "  where bedtools     # Windows",
      call. = FALSE
    )
  }

  if (!quiet) message("bedtools found at: ", bt_path)
  invisible(bt_path)
}
