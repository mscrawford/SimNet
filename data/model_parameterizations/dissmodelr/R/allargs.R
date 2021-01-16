#' Get function arguments with evaluated values
#'
#' @return List with function arguments
#' @export

allargs <- function() {
  # simplified from:
  # https://stackoverflow.com/questions/11885207/get-all-parameters-as-list/11892680
  # get formals for parent function
  parent_formals <- formals(sys.function(sys.parent(n = 1)))

  # Get names of implied arguments
  fnames <- names(parent_formals)

  # # Remove '...' from list of parameter names if it exists
  # fnames <- fnames[-which(fnames == '...')]

  # Get currently set values for named variables in the parent frame
  args <- evalq(as.list(environment()), envir = parent.frame())

  # # Get the list of variables defined in '...'
  # args <- c(args[fnames], evalq(list(...), envir = parent.frame()))

  args
}
