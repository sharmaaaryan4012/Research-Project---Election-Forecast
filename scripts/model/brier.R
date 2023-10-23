election_results <- list(
  AL = 0, AK = 0, AZ = 0, AR = 0, CA = 1,
  CO = 1, CT = 1, DE = 1, FL = 0, GA = 0,
  HI = 1, ID = 0, IL = 1, IN = 0, IA = 0,
  KS = 0, KY = 0, LA = 0, ME = 1, MD = 1,
  MA = 1, MI = 0, MN = 1, MS = 0, MO = 0,
  MT = 0, NE = 0, NV = 1, NH = 1, NJ = 1,
  NM = 1, NY = 1, NC = 0, ND = 0, OH = 0,
  OK = 0, OR = 1, PA = 0, RI = 1, SC = 0,
  SD = 0, TN = 0, TX = 0, UT = 0, VT = 1,
  VA = 1, WA = 1, WV = 0, WI = 0, WY = 0
)

predicted_probs <- c(0.974, 0.947, 0.905, 0.869, 0.861, 0.808, 0.753, 0.817, 0.739, 0.737, 0.649, 0.439, 0.223, 0.151, 0.0747, 0.0787, 0.0147, 0.00333, 0.004, 0.002)
states <- c("ME", "NM", "MN", "VA", "WI", "MI", "CO", "--", "PA", "NH", "NV", "FL", "OH", "IA", "AZ", "GA", "SC", "TX", "MO", "MS")
actual_outcomes <- unlist(election_results[states], use.names = FALSE)
brier_score <- mean((predicted_probs - actual_outcomes)^2)
print(brier_score)