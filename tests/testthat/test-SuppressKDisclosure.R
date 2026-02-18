
printInc <- FALSE
print_frames <- FALSE

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
  
  d$freq2 = round(d$freq/7)
  
  suppsums <- integer(0)
  
  sensitive <- vector("list", 2)
  sensitive[[2]] <- list(region = c("A", "C", "G"), main_income = c("pensions", "wages"))
  for (extend0 in c(TRUE, FALSE)) {
    for (i in 1:2) {
      for (singletonMethod in  c("anySumNOTprimary", "anySum0") ) {
        a <- SuppressKDisclosure(d, dimVar = 1:4, freqVar = "freq2", coalition = 3, 
                                 extend0 = extend0, sensitive = sensitive[[i]],
                                 whenEmptyUnsuppressed = NULL,
                                 singletonMethod = singletonMethod , 
                                 printInc = printInc)
        suppsums <- c(suppsums, sum(a$suppressed))
      }
    }
  }
  
  expect_identical(suppsums, c(63L, 69L, 38L, 38L, 50L, 52L, 33L, 41L))
  
  
  mm <- SSBtools::ModelMatrix(d, dimVar = 1:4, crossTable = TRUE)
  
  targ <- default_targeting(crossTable = mm$crossTable, x = mm$modelMatrix, 
                            sensitive = sensitive[[2]], 
                            identifying = list(region = "!8", main_income = "*"))
  
  targ[[1]][10:11, 1] <- "99"   # to test matching
  targ[[2]][10:11, 1] <- "999"  # 
  
  o1 <- SuppressKDisclosure(d, dimVar = 1:4, freqVar = "freq", 
                            coalition = 55, extend0 = FALSE, 
                            targeting = targ, whenEmptyUnsuppressed = NULL,
                            printInc = printInc, print_frames = print_frames, output = "all")
  
  me <- Matrix::Matrix(FALSE, nrow(targ$sensitive), nrow(targ$identifying))
  
  me[targ$sensitive$region == "A" & targ$sensitive$main_income == "assistance", 
     targ$identifying$region == "1" & targ$identifying$main_income ==  "Total"] <- TRUE
  
  
  targ$exclude_relations <- me
  
  
  o2 <- SuppressKDisclosure(d, dimVar = 1:4, freqVar = "freq", coalition = 55, 
                            extend0 = FALSE, targeting = targ, 
                            whenEmptyUnsuppressed = NULL,
                            printInc = printInc, print_frames = print_frames, 
                            output = "all")
  
  expect_identical(c(ncol(o2$xExtraPrimary), ncol(o1$xExtraPrimary)), 23:24)
  
  
  targ$targeting_exclude <- list(
    list(sensitive = list(region = "A", main_income = "assistance"),
         identifying = list(region = "1", main_income = "Total"))
  )
  
  o3 <- SuppressKDisclosure(d, dimVar = 1:4, freqVar = "freq", coalition = 55, 
                            extend0 = FALSE,
                            sensitive = sensitive[[2]], 
                            identifying = list(region = "!8", main_income = "*"),
                            whenEmptyUnsuppressed = NULL,
                            printInc = printInc, print_frames = print_frames, 
                            output = "all")
  
  identical(as.vector(table(SSBtools::DummyDuplicated(cbind(o1$xExtraPrimary, o3$xExtraPrimary), rnd = TRUE))),
            c(24L, 24L))
  
  
  o4 <- SuppressKDisclosure(d, dimVar = 1:4, freqVar = "freq", coalition = 55, 
                            extend0 = FALSE,
                            sensitive = sensitive[[2]], 
                            identifying = list(region = "!8", main_income = "*"),
                            targeting_exclude = list(
                              list(sensitive = list(region = "A", main_income = "assistance"),
                                   identifying = list(region = "1", main_income = "Total"))),
                            whenEmptyUnsuppressed = NULL,
                            printInc = printInc, print_frames = print_frames, 
                            output = "all")
  
  identical(as.vector(table(SSBtools::DummyDuplicated(cbind(o2$xExtraPrimary, o4$xExtraPrimary), rnd = TRUE))),
            c(23L, 22L))
  
  
  o5 <- SuppressKDisclosure(d, dimVar = 1:4, freqVar = "freq", coalition = 55, 
                            extend0 = FALSE,
                            sensitive = sensitive[[2]], 
                            identifying = list(region = "!8", main_income = "*"),
                            targeting_exclude = list(
                              list(sensitive = data.frame(region = "A", main_income = "assistance"),
                                   identifying = data.frame(region = "1", main_income = "Total"))),
                            whenEmptyUnsuppressed = NULL,
                            printInc = printInc, print_frames = print_frames, 
                            output = "all")
  
  identical(as.vector(table(SSBtools::DummyDuplicated(cbind(o2$xExtraPrimary, o5$xExtraPrimary), rnd = TRUE))),
            c(23L, 23L))
  
  
  
  
  ### Tests based on advanced examples using `targeting_exclude` and `targeting_include` 
  
  
  # Create a wrapper function to avoid repeating common arguments                                
  fun <- function(..., coalition = 7) {
    SuppressKDisclosure(SSBtoolsData("d3"), 
                        formula = ~(region + county)*main_income + region*months + county*main_income*months, 
                        freqVar = "freq", coalition = coalition ,
                        mc_hierarchies = list(main_income = c("special = assistance + other", 
                                                              "ordinary = pensions + wages")),
                        printInc = printInc, print_frames = print_frames, output = "all",
                        ...)}
  
  
  # Only the categories "assistance" and "wages" are considered sensitive
  # Also use "special" and "ordinary" as identifying categories (instead of "Total")
  a4 <- fun(sensitive = list(main_income = c("assistance", "wages")), 
            identifying = list(region = "*", main_income = c("special", "ordinary")))
  
  
  a4_ <- fun(targeting_include = list(list(
            sensitive = list(main_income = c("assistance", "wages")), 
            identifying = list(region = "*", main_income = c("special", "ordinary")))))
  
  expect_identical(a4, a4_)
  
  
  # As above, but additionally exclude regions i and j via the sensitive specification          
  a5 <- fun(sensitive = list(main_income = c("assistance", "wages")), 
            identifying = list(region = "*", main_income = c("special", "ordinary")), 
            targeting_exclude = list(list(sensitive = list(region = c("i", "j")))))
  
  # Same exclusion as above, but specified via identifying instead of sensitive
  # Here `main_income` must also be specified, since the default for identifying is "Total" 
  a6 <- fun(sensitive = list(main_income = c("assistance", "wages")), 
            identifying = list(region = "*", main_income = c("special", "ordinary")), 
            targeting_exclude = list(list(identifying = list(region = c("i", "j"), 
                                                             main_income = "*"))))
  
  
  
  # As above, but use a data.frame for precise specification of relations
  # Therefore, "V ordinary–pensions" is no longer included                                     
  a8 <- fun(sensitive = list(main_income = c("assistance", "wages")), 
            identifying = list(region = "*", main_income = c("special", "ordinary")), 
            targeting_exclude = list(list(identifying = list(region = c("i", "j"), 
                                                             main_income = "*"))), 
            targeting_include = list(
              list(identifying = data.frame(region = c("14", "U", "V", "X"), 
                                            main_income = c("special", "ordinary"), 
                                            months = c("m10m12", "Total")), 
                   sensitive = list(region = c("m01m05"), 
                                    main_income = c("pensions", "assistance")))))    
  
  # Specify the same relations as above, but in a different way
  # Using multiple list elements                                    
  a9 <- fun(sensitive = list(main_income = c("assistance", "wages")), 
            identifying = list(region = "*", main_income = c("special", "ordinary")), 
            targeting_exclude = list(list(identifying = list(region = c("i", "j"), 
                                                             main_income = "*"))), 
            targeting_include = list(
              list(identifying = list(region = "14", 
                                      main_income = "special", 
                                      months = "m10m12"), 
                   sensitive = list(region = "14", 
                                    main_income = "assistance", 
                                    months = "m10m12")), 
              list(identifying = list(region = c("U", "X"), 
                                      main_income = "ordinary", 
                                      months = "Total"), 
                   sensitive = list(region = c("U", "X"), 
                                    main_income = "pensions", 
                                    months = "Total")))) 
  
  expect_identical(c(ncol(a5$xExtraPrimary), ncol(a6$xExtraPrimary), 
                     ncol(a8$xExtraPrimary), ncol(a9$xExtraPrimary)),
                   c(22L, 22L, 25L, 25L))
  
  expect_identical(a5$publish, a6$publish)
  expect_identical(a8$publish, a9$publish)
  
})
