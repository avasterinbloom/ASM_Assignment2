#Problem 1

df_eq <- read.csv("https://www.isibang.ac.in/~rsen/ASM1/EarthquakeData.csv")

time_to_fraction <- function(t) {
  if(is.na(t)) return(NA)
  t_int <- floor(t)
  dec <- t - t_int
  ss <- t_int %% 100 + dec
  mm <- (t_int %/% 100) %% 100
  hh <- (t_int %/% 10000)
  return((hh * 3600 + mm * 60 + ss) / 86400)
}


df_eq$frac_day <- sapply(df_eq$Time.hhmmss.mm.UTC., time_to_fraction)
t_data <- sort(na.omit(df_eq$frac_day))
n <- length(t_data)

# 1(a) Plotting Empirical CDF vs Uniform CDF
plot(ecdf(t_data), main="Empirical vs Uniform CDF", 
     xlab="Fraction of Day (t)", ylab="CDF", col="blue")
curve(punif(x), add=TRUE, col="red", lty=2)
legend("bottomright", legend=c("Empirical CDF", "Uniform CDF"), 
       col=c("blue", "red"), lty=c(1, 2))

# 1(b) Plot of Gn(t)
F_hat <- ecdf(t_data)(t_data)
Gn <- sqrt(n) * (F_hat - t_data)
plot(t_data, Gn, type="l", col="purple", main="Plot of G_n(t)", 
     xlab="Fraction of Day (t)", ylab="G_n(t)")
abline(h=0, col="black", lty=2)

# 1(c) Kolmogorov-Smirnov Test
ks_result <- ks.test(t_data, "punif")
print(ks_result)



#Problem 2
germ_prod <- read.table(https://www.isibang.ac.in/~rsen/ASM1/GermProd.txt, header = TRUE)

nw_kernel <- function(x_eval, X, Y, h, kernel_type="epanechnikov") {
  preds <- numeric(length(x_eval))
  for (i in seq_along(x_eval)) {
    u <- (X - x_eval[i]) / h
    if (kernel_type == "gaussian") {
      K <- dnorm(u)
    } else if (kernel_type == "epanechnikov") {
      K <- 0.75 * (1 - u^2) * (abs(u) <= 1)
    } else if (kernel_type == "uniform") {
      K <- 0.5 * (abs(u) <= 1)
    }
    preds[i] <- sum(K * Y) / sum(K)
  }
  return(preds)
}

x_grid <- seq(min(germ_prod$month), max(germ_prod$month), length.out=500)

# 2(a) Bandwidth 24 with Gaussian, Epanechnikov, and Uniform kernels
y_gauss <- nw_kernel(x_grid, germ_prod$month, germ_prod$production, 24, "gaussian")
y_epan  <- nw_kernel(x_grid, germ_prod$month, germ_prod$production, 24, "epanechnikov")
y_unif  <- nw_kernel(x_grid, germ_prod$month, germ_prod$production, 24, "uniform")

plot(germ_prod$month, germ_prod$production, col="gray", 
     xlab="Month", ylab="Production", main="Kernel Regression (h=24)")
lines(x_grid, y_gauss, col="red", lwd=2)
lines(x_grid, y_epan, col="blue", lwd=2)
lines(x_grid, y_unif, col="green", lwd=2)
legend("topleft", legend=c("Gaussian", "Epanechnikov", "Uniform"), 
       col=c("red", "blue", "green"), lwd=2)

# 2(b) Epanechnikov Kernel with Bandwidths 10, 20, 30
y_bw10 <- nw_kernel(x_grid, germ_prod$month, germ_prod$production, 10, "epanechnikov")
y_bw20 <- nw_kernel(x_grid, germ_prod$month, germ_prod$production, 20, "epanechnikov")
y_bw30 <- nw_kernel(x_grid, germ_prod$month, germ_prod$production, 30, "epanechnikov")

plot(germ_prod$month, germ_prod$production, col="gray", 
     xlab="Month", ylab="Production", main="Epanechnikov Kernel (Varying Bandwidths)")
lines(x_grid, y_bw10, col="red", lwd=2)
lines(x_grid, y_bw20, col="blue", lwd=2)
lines(x_grid, y_bw30, col="green", lwd=2)
legend("topleft", legend=c("h=10", "h=20", "h=30"), 
       col=c("red", "blue", "green"), lwd=2)
