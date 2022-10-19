
data("RiskManagerCost")

expect_error(robustbetareg(SIZELOG ~ FIRMCOST + INDCOST,
                           data = RiskManagerCost, type = "LMDPDE", alpha = 0,
                           link.phi = "log"))

expect_error(robustbetareg(SIZELOG ~ FIRMCOST + INDCOST,
                           data = RiskManagerCost, type = "LMDPDE", alpha = -3,
                           link.phi = "log"))

expect_warning(dEGB(0.2, mu = -0.3, phi = 1))

