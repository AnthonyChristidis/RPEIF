# -----------------------------------------
# Test Script - Warnings for returns data
# -----------------------------------------
                  
# Required libraries
library(IFs)

# Context of test script
context("Verify wrapper function.")

# The wrapper function and the risk function should return the same output for an IF TS
test_that("Equal output for IF TS using wrapper function.", {
  data(edhec, package="PerformanceAnalytics")
  colnames(edhec) = c("CA", "CTAG", "DIS", "EM","EMN", "ED", "FIA",
                      "GM", "LS", "MA", "RV", "SS", "FoF") 
  expect_equal(IF(risk="mean", returns=edhec[,"CA"]), IF.mean(returns=edhec[,"CA"]))
  expect_equal(IF(risk="SD", returns=edhec[,"CA"]), IF.SD(returns=edhec[,"CA"]))
  expect_equal(IF(risk="SSD", returns=edhec[,"CA"]), IF.SSD(returns=edhec[,"CA"]))
  expect_equal(IF(risk="LPM", returns=edhec[,"CA"]), IF.LPM(returns=edhec[,"CA"]))
  expect_equal(IF(risk="ES", returns=edhec[,"CA"]), IF.ES(returns=edhec[,"CA"]))
  expect_equal(IF(risk="VaR", returns=edhec[,"CA"]), IF.VaR(returns=edhec[,"CA"]))
  expect_equal(IF(risk="SR", returns=edhec[,"CA"]), IF.SR(returns=edhec[,"CA"]))
  expect_equal(IF(risk="SoR", returns=edhec[,"CA"]), IF.SoR(returns=edhec[,"CA"]))
  expect_equal(IF(risk="ESratio", returns=edhec[,"CA"]), IF.ESratio(returns=edhec[,"CA"]))
  expect_equal(IF(risk="VaRratio", returns=edhec[,"CA"]), IF.VaRratio(returns=edhec[,"CA"]))
  expect_equal(IF(risk="RachR", returns=edhec[,"CA"]), IF.RachR(returns=edhec[,"CA"]))
  expect_equal(IF(risk="Omega", returns=edhec[,"CA"]), IF.Omega(returns=edhec[,"CA"]))
})
