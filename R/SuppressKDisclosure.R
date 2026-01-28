#' K-disclosure suppression
#'
#' A function for suppressing frequency tables using the k-disclosure method.
#' 
#' @details
#' The `targeting` specification is a named list that may contain the following
#' optional elements. References to `crossTable` below refer to a data frame
#' that may be extended after applying `mc_hierarchies`.
#'
#' \describe{
#'   \item{`identifying`}{A data frame containing selected rows from
#'   `crossTable`. Membership in the cells represented by these rows is regarded
#'   as information that an intruder may already know. If omitted, it defaults
#'   to `crossTable`.}
#'
#'   \item{`disclosive`}{A data frame containing selected rows from
#'   `crossTable`. If an intruder can infer membership in the cells represented
#'   by these rows, this is considered an unacceptable disclosure, subject to
#'   any further specification provided by `is_disclosive`. If omitted, it
#'   defaults to `crossTable`.}
#'
#'   \item{`is_disclosive`}{A data frame with the same structure as
#'   `disclosive`, but with logical variables. It indicates which codes in
#'   `disclosive` are regarded as disclosive. When specified, disclosure is
#'   assessed by which codes within a revealed cell are disclosed. If omitted,
#'   it is equivalent to a data frame where all elements are `TRUE`.}
#' }
#'
#' The argument `targeting` may also be a function that returns such a list.
#' This works similarly to supplied functions in `GaussSuppressionFromData()`.
#' Note, however, that the function operates on possibly extended versions of
#' `freq`, `x`, and `crossTable` that reflect the use of `mc_hierarchies`, when
#' applicable.
#'
#' @param data a data.frame representing the data set
#' @param coalition numeric vector of length one, representing possible size of an
#' attacking coalition. This parameter corresponds to the parameter k in the
#' definition of k-disclosure.
#' @param dimVar The main dimensional variables and additional aggregating
#' variables. This parameter can be  useful when hierarchies and formula are
#' unspecified.
#' @param formula A model formula
#' @param hierarchies List of hierarchies, which can be converted by
#' \code{\link[SSBtools]{AutoHierarchies}}. Thus, the variables can also be coded by
#' `"rowFactor"` or `""`, which correspond to using the categories in the data.
#' @param freqVar name of the frequency variable in `data`
#' @param mc_hierarchies a hierarchy representing meaningful combinations to be
#' protected. Default value is `NULL`.
#' @param upper_bound numeric value representing minimum count considered safe.
#' Default set to `Inf`
#' @param ... parameters passed to children functions
#' @inheritParams GaussSuppressionFromData
#'
#' @param targeting NULL, a list, or a function that returns a list specifying
#' attribute disclosure scenarios. See Details.
#' 
#' @param print_frames Logical. If TRUE, additional data frames are printed to
#' the console. When `mc_hierarchies` is used, this includes a
#' data frame with hidden results. In addition, a data frame containing the
#' primary suppressed difference cells is printed. The default is FALSE.

#'
#' @return A data.frame containing the publishable data set, with a boolean
#' variable `$suppressed` representing cell suppressions.
#' @export
#'
#' @author Daniel P. Lupp and Øyvind Langsrud
#'
#' @examples
#' # data
#' data <- SSBtools::SSBtoolsData("mun_accidents")
#'
#' # hierarchies as DimLists
#' mun <- data.frame(levels = c("@@", rep("@@@@", 6)),
#' codes = c("Total", paste("k", 1:6, sep = "")))
#' inj <- data.frame(levels = c("@@", "@@@@" ,"@@@@", "@@@@", "@@@@"),
#' codes = c("Total", "serious", "light", "none", "unknown"))
#' dimlists <- list(mun = mun, inj = inj)
#'
#' inj2 <- data.frame(levels = c("@@", "@@@@", "@@@@@@" ,"@@@@@@", "@@@@", "@@@@"),
#' codes = c("Total", "injured", "serious", "light", "none", "unknown"))
#' inj3 <- data.frame(levels = c("@@", "@@@@", "@@@@" ,"@@@@", "@@@@"),
#' codes = c( "shadowtotal", "serious", "light", "none", "unknown"))
#' mc_dimlist <- list(inj = inj2)
#' mc_nomargs <- list(inj = inj3)
#'
#' #' # Example with formula, no meaningful combination
#' out <- SuppressKDisclosure(data, coalition = 1, freqVar = "freq", formula = ~mun*inj)
#'
#' # Example with hierarchy and meaningful combination
#' out2 <- SuppressKDisclosure(data, coalition = 1, freqVar = "freq",
#' hierarchies = dimlists, mc_hierarchies = mc_dimlist)
#'
#' #' # Example of table without mariginals, and mc_hierarchies to protect
#' out3 <- SuppressKDisclosure(data, coalition = 1, freqVar = "freq",
#' formula = ~mun:inj, mc_hierarchies = mc_nomargs )
SuppressKDisclosure <- function(data,
                                coalition = 0,
                                mc_hierarchies = NULL,
                                upper_bound = Inf,
                                dimVar = NULL,
                                formula = NULL,
                                hierarchies = NULL,
                                freqVar = NULL,
                                targeting = NULL,
                                print_frames = FALSE,
                                ...,
                                spec = PackageSpecs("kDisclosureSpec")) {
  additional_params <- list(...)
  if (length(additional_params)) {
    if ("singletonMethod" %in% names(additional_params) &
        "none" %in% additional_params[["singletonMethod"]])
      warning(
        "SuppressKDisclosure should use a singleton method for protecting the zero singleton problem. The output might not be safe, consider rerunning with a singleton method (default)."
      )
  }
  GaussSuppressionFromData(
    data,
    hierarchies = hierarchies,
    formula = formula,
    dimVar = dimVar,
    freqVar = freqVar,
    coalition = coalition,
    mc_hierarchies = mc_hierarchies,
    upper_bound = upper_bound,
    spec = spec,
    targeting = targeting,
    print_frames = print_frames,
    ...
  )
}

#' Construct primary suppressed difference matrix
#'
#' Function for constructing model matrix columns representing primary suppressed
#' difference cells
#'
#' @inheritParams SuppressKDisclosure
#' @param crossTable crossTable generated by parent function
#' @param x ModelMatrix generated by parent function
#'
#' @return dgCMatrix corresponding to primary suppressed cells
#' @export
#'
#' @author Daniel P. Lupp and Øyvind Langsrud
KDisclosurePrimary <- function(data,
                               x,
                               crossTable,
                               freqVar,
                               mc_hierarchies = NULL,
                               coalition = 1,
                               upper_bound = Inf,
                               targeting = NULL,
                               print_frames = FALSE,
                               ...) {
  
  
  
  mc_obj <- X_from_mc(
    data = data,
    x = x,
    crossTable = crossTable,
    mc_hierarchies = mc_hierarchies,
    freqVar = freqVar,
    coalition = coalition,
    upper_bound = upper_bound,
    returnNewCrossTable = TRUE,
    ...
  )
  
  orig_nrow_crossTable <- nrow(crossTable)
  
  x <- cbind(x, mc_obj$x)
  crossTable <- rbind(crossTable, mc_obj$crossTable)
  
  freq <- as.vector(crossprod(x, data[[freqVar]]))
  
  
  if (print_frames & !is.null(mc_obj)) {
    r <- SSBtools::SeqInc(orig_nrow_crossTable + 1, nrow(crossTable))
    hidden_cells <- cbind(crossTable[r, ,drop = FALSE], freq = freq[r])
    rownames(hidden_cells) <- NULL
    cat("\n----- hidden cells from mc_hierarchies -----\n")
    print(hidden_cells)
  }
  
  
  if(is.function(targeting)) {
    targeting <- targeting(..., freq = freq, x = x, crossTable = crossTable)
  }
  
  identifying <- targeting$identifying
  disclosive  <- targeting$disclosive
  is_disclosive <- targeting$is_disclosive
  
  use_is_disclosive <- !is.null(is_disclosive)
  
  if (use_is_disclosive | print_frames) {
    if (is.null(identifying)) {
      identifying <- crossTable
    }
    if (is.null(disclosive)) {
      disclosive <- crossTable
    }
  }
    
  if (use_is_disclosive) {
    validate_is_disclosive(is_disclosive, disclosive)
    
    if (isTRUE(all(is_disclosive))) {
      is_disclosive <- NULL
      use_is_disclosive <- FALSE
    }
  }
  
  if (!is.null(identifying) | !is.null(disclosive)) {
    
    if (!is.null(identifying)) {
      ma <- SSBtools::Match(identifying, crossTable)
      ma <- ma[!is.na(ma)]
      if (use_is_disclosive | print_frames) {
        identifying <- identifying[!is.na(ma), ]
      }
      y <- x[, ma]
    } else {
      y <- x
    }
    if (!use_is_disclosive) {
      y <- y[, !SSBtools::DummyDuplicated(y, rnd = TRUE), drop = FALSE]
    }
    
    if (!is.null(disclosive)) {
      ma <- SSBtools::Match(disclosive, crossTable)
      ma <- ma[!is.na(ma)]
      if (use_is_disclosive| print_frames) {
        disclosive <- disclosive[!is.na(ma), ]
      }
      if (use_is_disclosive) {
        is_disclosive <- is_disclosive[!is.na(ma), ]
      }
      x <- x[, ma]
    }
    if (!use_is_disclosive) {
      x <- x[, !SSBtools::DummyDuplicated(x, rnd = TRUE), drop = FALSE]
    }
  } else {
    x <- x[, !SSBtools::DummyDuplicated(x, rnd = TRUE), drop = FALSE]
    y <- x
  }
  
  if (use_is_disclosive) {  # Extra check after modifications  
    validate_is_disclosive(is_disclosive, disclosive)
  }
  
  FindDifferenceCells(
    x = x,
    y = y,
    freq_x = as.vector(crossprod(x, data[[freqVar]])),
    freq_y = as.vector(crossprod(y, data[[freqVar]])),
    coalition = coalition,
    upper_bound = upper_bound,
    crossTable = crossTable,
    identifying = identifying,
    disclosive = disclosive,
    is_disclosive = is_disclosive,
    print_frames 
  )
}



FindDifferenceCells <- function(x,
                                y = x,
                                freq_x,
                                freq_y = freq_x,
                                coalition,
                                upper_bound = Inf,
                                crossTable, # used only for via nrow()
                                identifying,
                                disclosive,
                                is_disclosive,
                                print_frames = FALSE
                                ) {
  xty <- As_TsparseMatrix(crossprod(x, y))
  colSums_y_xty_j_1 <- colSums(y)[xty@j + 1]
  # finds children in x and parents in y
  r <- colSums(x)[xty@i + 1] == xty@x & 
    colSums_y_xty_j_1     != xty@x & 
    (colSums_y_xty_j_1 - xty@x) <= upper_bound
  
  if (!any(r)) {
    return(rep(FALSE, nrow(crossTable)))
  }
  
  child <- xty@i[r] + 1L
  parent <- xty@j[r] + 1L
  
  freq_diff <- freq_y[parent] - freq_x[child]
  
  disclosures <- 
    freq_x[child] > 0  &
    freq_y[parent] > 0  &
    freq_diff <= coalition
  
  if (!any(disclosures)) {
    return(rep(FALSE, nrow(crossTable)))
  }
  
  freq_diff <- freq_diff[disclosures]
  parent <- parent[disclosures]
  child <- child[disclosures]
  
  use_is_disclosive <- !is.null(is_disclosive)
  
  if (use_is_disclosive | print_frames) {
    identifying <- identifying[parent, , drop = FALSE]
    disclosive <- disclosive[child, , drop = FALSE]
  }
  
  if (use_is_disclosive) {
    is_disclosive <- as.matrix(is_disclosive)
    is_disclosive <- is_disclosive[child, , drop = FALSE]
    
    same_codes <- identifying == disclosive
    same_codes[!is_disclosive] <- TRUE
    same_row   <- apply(same_codes, 1, all)
    
    parent <- parent[!same_row]
    child <- child[!same_row]
    if (print_frames) {
      identifying <- identifying[!same_row, , drop = FALSE]
      disclosive <- disclosive[!same_row, , drop = FALSE]
      freq_diff <- freq_diff[!same_row]
    }
  }
  
  if (print_frames) {
    cat("\n---- primary suppressed difference cells ---\n")
    print_difference_cells(identifying, disclosive, freq_diff)
    cat("\n")
  }
  
  diff_matrix <- drop0(y[, parent, drop = FALSE] - 
                       x[, child, drop = FALSE])
  
  
  diff_matrix[, !SSBtools::DummyDuplicated(diff_matrix, rnd = TRUE), drop = FALSE]
  
}







# Written by ChatGPT
validate_is_disclosive <- function(is_disclosive, disclosive) {
  if (!is.data.frame(is_disclosive)) {
    stop("`is_disclosive` must be a data frame.", call. = FALSE)
  }
  
  if (!identical(dim(is_disclosive), dim(disclosive))) {
    stop(
      "`is_disclosive` must have the same dimensions as `disclosive`.",
      call. = FALSE
    )
  }
  
  if (!identical(names(is_disclosive), names(disclosive))) {
    stop(
      "`is_disclosive` must have the same variable names as `disclosive`.",
      call. = FALSE
    )
  }
  
  if (!all(vapply(is_disclosive, is.logical, logical(1)))) {
    stop(
      "All variables in `is_disclosive` must be logical.",
      call. = FALSE
    )
  }
  
  # no return value needed
}



print_difference_cells <- function(identifying, disclosive, freq_diff) {
  r <- identifying != disclosive
  identifying[r] <- paste(identifying[r], disclosive[r], sep = "-")
  identifying$diff <- freq_diff
  rownames(identifying) <- NULL
  print(identifying)
}
