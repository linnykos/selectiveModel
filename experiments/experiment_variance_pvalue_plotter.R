png("../figures/pvalue_stability.png", heigh = 2000, width = 2500, units = "px",
    res = 300)
par(mfrow = c(2,2))
hist(res_vec_known_0, breaks = seq(0, 1, length.out = 50), col = "gray",
     xlab = "P-value", main = "Delta = 0 (Null), Known sigma")
hist(res_vec_unknown_0, breaks = seq(0, 1, length.out = 50), col = "gray",
     xlab = "P-value", main = "Delta = 0 (Null), Unknown sigma")
hist(res_vec_known_signal, breaks = seq(0, 1, length.out = 50), col = "gray",
     xlab = "P-value", main = "Delta = 1 (Signal), Known sigma")
graphics.off()
