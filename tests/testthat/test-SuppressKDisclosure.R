
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
  
})
