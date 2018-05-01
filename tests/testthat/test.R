Sys.setenv("R_TESTS" = "")
library(testthat)
library(eggCountsExtra)

data(epgs)
t1 <- w_quant(epgs$after)
expect_that(t1, is_a("numeric"))

t2 <- w_ratio(epgs$before, epgs$after)
expect_that(t2, is_a("numeric"))

t3 <- getCode("Po")
expect_that(t3, is_a("character"))
