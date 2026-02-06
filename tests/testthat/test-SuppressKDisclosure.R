
printInc <- FALSE

options(GaussSuppression.action_unused_dots = "abort")

test_that("SuppressKDisclosure", {
  
  mun_accidents <- SSBtoolsData("mun_accidents")

  
  # hierarchies as DimLists
  mun <- data.frame(levels = c("@", rep("@@", 6)),
                    codes = c("Total", paste("k", 1:6, sep = "")))
  inj <- data.frame(levels = c("@", "@@" ,"@@", "@@", "@@"),
                    codes = c("Total", "serious", "light", "none", "unknown"))
  dimlists <- list(mun = mun, inj = inj)
  
  inj2 <- data.frame(levels = c("@", "@@", "@@@" ,"@@@", "@@", "@@"),
                     codes = c("Total", "injured", "serious", "light", "none", "unknown"))
  inj3 <- data.frame(levels = c("@", "@@", "@@" ,"@@", "@@"),
                     codes = c( "shadowtotal", "serious", "light", "none", "unknown"))
  mc_dimlist <- list(inj = inj2)
  mc_nomargs <- list(inj = inj3)
  
  #' # Example with formula, no meaningful combination
  out1 <- SuppressKDisclosure(mun_accidents, coalition = 1, freqVar = "freq", formula = ~mun*inj,
                              printInc = printInc)
  
  # Example with hierarchy and meaningful combination
  out2 <- SuppressKDisclosure(mun_accidents, coalition = 1, freqVar = "freq",
                              hierarchies = dimlists, mc_hierarchies = mc_dimlist,
                              printInc = printInc)
  
  #' # Example of table without mariginals, and mc_hierarchies to protect
  out3 <- SuppressKDisclosure(mun_accidents, coalition = 1, freqVar = "freq",
                              formula = ~mun:inj, mc_hierarchies = mc_nomargs,
                              printInc = printInc)
  
  expect_identical(as.list(table(out1[out1[["suppressed"]], "inj"])), 
                   list(light = 2L, none = 4L, serious = 4L))
  expect_identical(as.list(table(out2[out2[["suppressed"]], "inj"])), 
                   list(light = 7L, none = 2L, serious = 4L, unknown = 5L))
  expect_identical(as.list(table(out3[out3[["suppressed"]], "inj"])), 
                   list(none = 1L, unknown = 3L))
  
  
  
  d2 <- SSBtoolsData("d2")
  
  minus <- c(5, 7, 9, 11, 12, 13, 21, 26, 28, 29, 31, 34, 37, 38)
  d <- d2[-minus, ]
  
  suppsums <- integer(0)
  
  disclosive <- vector("list", 2)
  disclosive[[2]] <- list(region = c("A", "C", "G"), main_income = c("pensions", "wages"))
  for (extend0 in c(TRUE, FALSE)) {
    for (i in 1:2) {
      a <- SuppressKDisclosure(d, dimVar = 1:4, freqVar = "freq", coalition = 3, 
                               extend0 = extend0, disclosive = disclosive[[i]],
                               whenEmptyUnsuppressed = NULL,
                               printInc = printInc)
      suppsums <- c(suppsums, sum(a$suppressed))
    }
  }
  expect_identical(suppsums, c(53L, 28L, 32L, 22L))
  
})
