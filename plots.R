#Script for generating plots in Figure 1 of Cattaneo, Chandak, Jansson and Ma (2022)

#install companion package
install.packages("lpcde")

#load companion package
library("lpcde")

# set seed for replicability
set.seed(58)

# setting parameters
n = 1000
ng = 20
y_grid = seq(-2, 2, length.out=ng)
x = 0
CIsimul = 2000
alpha = 0.05

#CDF PLOTS
order = "cdf"
type = "wbc"

mu = 0
p = 2
true_dens = pnorm(y_grid)

# initializing
est = matrix(0L, nrow=ng, ncol=1)
ci_l = matrix(0L, nrow=ng, ncol=1)
ci_u = cb_l = cb_u = ci_l

# simulate data
x_data = as.matrix(rnorm(5*n))
x_data = x_data[abs(x_data)<=3.2][1:n]
y_data = as.matrix(rnorm(5*n))
y_data = y_data[abs(y_data)<=3.2][1:n]

# estimating bandwidth
bw = lpcde::lpbwcde(y_data=y_data, x_data=x_data, x=x, y_grid = y_grid,
                    p=p, mu=mu, kernel_type="epanechnikov", bw_type="imse-rot")$BW[,2]

# estimating density and confidence bands
model1 = lpcde::lpcde(y_data=y_data, x_data=x_data, x=x, y_grid = y_grid,
                      mu=mu, p=p, bw=bw, kernel_type = "epanechnikov", var_type="ustat")

est = est + model1$Estimate[,3]

# confidence bands
corrMat = sqrt(abs(sweep(sweep(model1$CovMat$CovMat, MARGIN=1, FUN="*",
                               STATS=1/model1$Estimate[, "se"]), MARGIN=2, FUN="*",
                         STATS=1/model1$Estimate[, "se"])))
normalSimu = try(
     MASS::mvrnorm(n=CIsimul, mu=rep(0,nrow(corrMat)),
                   Sigma=Matrix::nearPD(corrMat)$mat), silent=TRUE)
z_val = stats::quantile(apply(normalSimu, MARGIN=1,
                              FUN=function(x) {max(abs(x))}), 1 - alpha)
ci_l = ci_l + (model1$Estimate[,3] - 1.96*sqrt(abs(model1$Estimate[,5])))
ci_u = ci_u +(model1$Estimate[,3] + 1.96*sqrt(abs(model1$Estimate[,5])))
cb_l = cb_l + (model1$Estimate[,3] - z_val*sqrt(abs(model1$Estimate[,5])))
cb_u = cb_u + (model1$Estimate[,3] + z_val*sqrt(abs(model1$Estimate[,5])))

# keeping only alternate confidnece intervals
for (i in 1:ng){
  if (i%%2==0){
    ci_l[i]=NA
    ci_u[i]=NA
  }
}

# creating dataframe
dataset = data.frame(cbind(y_grid, est, true_dens, ci_l, ci_u, cb_l, cb_u))
colnames(dataset) = c("y_grid", "est", "true_dens", "ci_l", "ci_u", "cb_l", "cb_u")

# plotting
temp_plot = ggplot2::ggplot() + ggplot2::theme_bw()
temp_plot = temp_plot+ ggplot2::ylim(-0.5, 1.4)
temp_plot = temp_plot + ggplot2::labs(x="", y="CDF") + ggplot2::ggtitle("")
temp_plot = temp_plot +
  ggplot2::geom_errorbar(data=dataset, ggplot2::aes(x=y_grid, ymin=ci_l, ymax=ci_u),
                         alpha=1, col=1, width=0.1)
temp_plot = temp_plot + ggplot2::geom_line(data=dataset,
                                           ggplot2::aes(x=y_grid, y=est), lwd=0.8)
temp_plot = temp_plot +
  ggplot2::geom_line(data=dataset, ggplot2::aes(x=y_grid, y=true_dens), col=2, lwd=0.8)
temp_plot = temp_plot +
  ggplot2::geom_ribbon(data=dataset, ggplot2::aes(x=y_grid, ymin=cb_l, ymax=cb_u),
                       alpha=0.2, fill=1)
temp_plot


type = "rbc"
mu = 0
p = 2
true_dens = pnorm(y_grid)

# initializing
est = matrix(0L, nrow=ng, ncol=1)
ci_l = matrix(0L, nrow=ng, ncol=1)
ci_u = cb_l = cb_u = ci_l

# simulate data
x_data = as.matrix(rnorm(5*n))
x_data = x_data[abs(x_data)<=3.2][1:n]
y_data = as.matrix(rnorm(5*n))
y_data = y_data[abs(y_data)<=3.2][1:n]

# estimating bandwidth
bw = lpcde::lpbwcde(y_data=y_data, x_data=x_data, x=x, y_grid = y_grid,
                    p=p, mu=mu, kernel_type="epanechnikov", bw_type="imse-rot")$BW[,2]

# estimating density and confidence bands
model1 = lpcde::lpcde(y_data=y_data, x_data=x_data, x=x, y_grid = y_grid,
                      mu=mu, p=p, bw=bw, kernel_type = "epanechnikov", var_type="ustat")

est = est + model1$Estimate[,3]

# confidence bands
corrMat = sqrt(abs(sweep(sweep(model1$CovMat$CovMat_RBC, MARGIN=1, FUN="*",
                               STATS=1/model1$Estimate[, "se_RBC"]), MARGIN=2, FUN="*",
                         STATS=1/model1$Estimate[, "se_RBC"])))
normalSimu = try(
  MASS::mvrnorm(n=CIsimul, mu=rep(0,nrow(corrMat)),
                Sigma=Matrix::nearPD(corrMat)$mat), silent=TRUE)
z_val = stats::quantile(apply(normalSimu, MARGIN=1,
                              FUN=function(x) {max(abs(x))}), 1 - alpha)
print(z_val)
ci_l = ci_l + (model1$Estimate[,4] - 1.96*sqrt(abs(model1$Estimate[,6])))
ci_u = ci_u + (model1$Estimate[,4] + 1.96*sqrt(abs(model1$Estimate[,6])))
cb_l = cb_l + (model1$Estimate[,4] - z_val*sqrt(abs(model1$Estimate[,6])))
cb_u = cb_u + (model1$Estimate[,4] + z_val*sqrt(abs(model1$Estimate[,6])))

# keeping only alternate confidnece intervals
for (i in 1:ng){
  if (i%%2==0){
    ci_l[i]=NA
    ci_u[i]=NA
  }
}

# creating dataframe
dataset = data.frame(cbind(y_grid, est, true_dens, ci_l, ci_u, cb_l, cb_u))
colnames(dataset) = c("y_grid", "est", "true_dens", "ci_l", "ci_u", "cb_l", "cb_u")

# plotting
temp_plot = ggplot2::ggplot() + ggplot2::theme_bw()
temp_plot = temp_plot+ ggplot2::ylim(-0.5, 1.4)
temp_plot = temp_plot + ggplot2::labs(x="", y="CDF") + ggplot2::ggtitle("")
temp_plot = temp_plot +
  ggplot2::geom_errorbar(data=dataset, ggplot2::aes(x=y_grid, ymin=ci_l, ymax=ci_u),
                         alpha=1, col=1, width=0.1)
temp_plot = temp_plot + ggplot2::geom_line(data=dataset,
                                           ggplot2::aes(x=y_grid, y=est), lwd=0.8)
temp_plot = temp_plot +
  ggplot2::geom_line(data=dataset, ggplot2::aes(x=y_grid, y=true_dens), col=2, lwd=0.8)
temp_plot = temp_plot +
  ggplot2::geom_ribbon(data=dataset, ggplot2::aes(x=y_grid, ymin=cb_l, ymax=cb_u),
                       alpha=0.2, fill=1)
temp_plot



#PDF PLOTS
type = "wbc"
mu = 1
p = 2
true_dens = dnorm(y_grid)

# initializing
est = matrix(0L, nrow=ng, ncol=1)
ci_l = matrix(0L, nrow=ng, ncol=1)
ci_u = cb_l = cb_u = ci_l

# simulate data
x_data = as.matrix(rnorm(5*n))
x_data = x_data[abs(x_data)<=3.2][1:n]
y_data = as.matrix(rnorm(5*n))
y_data = y_data[abs(y_data)<=3.2][1:n]

# estimating bandwidth
bw = lpcde::lpbwcde(y_data=y_data, x_data=x_data, x=x, y_grid = y_grid,
                    p=p, mu=mu, kernel_type="epanechnikov", bw_type="imse-rot")$BW[,2]

# estimating density and confidence bands
model1 = lpcde::lpcde(y_data=y_data, x_data=x_data, x=x, y_grid = y_grid,
                      mu=mu, p=p, bw=bw, kernel_type = "epanechnikov", var_type="ustat")

est = est + model1$Estimate[,3]

# confidence bands
corrMat = sqrt(abs(sweep(sweep(model1$CovMat$CovMat, MARGIN=1, FUN="*",
                               STATS=1/model1$Estimate[, "se"]), MARGIN=2, FUN="*",
                         STATS=1/model1$Estimate[, "se"])))
normalSimu = try(
     MASS::mvrnorm(n=CIsimul, mu=rep(0,nrow(corrMat)),
                   Sigma=Matrix::nearPD(corrMat)$mat), silent=TRUE)
z_val = stats::quantile(apply(normalSimu, MARGIN=1,
                              FUN=function(x) {max(abs(x))}), 1 - alpha)
ci_l = ci_l + (model1$Estimate[,3] - 1.96*sqrt(abs(model1$Estimate[,5])))
ci_u = ci_u +(model1$Estimate[,3] + 1.96*sqrt(abs(model1$Estimate[,5])))
cb_l = cb_l + (model1$Estimate[,3] - z_val*sqrt(abs(model1$Estimate[,5])))
cb_u = cb_u + (model1$Estimate[,3] + z_val*sqrt(abs(model1$Estimate[,5])))

# keeping only alternate confidnece intervals
for (i in 1:ng){
  if (i%%2==0){
    ci_l[i]=NA
    ci_u[i]=NA
  }
}

# creating dataframe
dataset = data.frame(cbind(y_grid, est, true_dens, ci_l, ci_u, cb_l, cb_u))
colnames(dataset) = c("y_grid", "est", "true_dens", "ci_l", "ci_u", "cb_l", "cb_u")

# plotting
temp_plot = ggplot2::ggplot() + ggplot2::theme_bw()
temp_plot = temp_plot+ ggplot2::ylim(-0.5, 0.7)
temp_plot = temp_plot + ggplot2::labs(x="", y="Density") + ggplot2::ggtitle("")
temp_plot = temp_plot +
  ggplot2::geom_errorbar(data=dataset, ggplot2::aes(x=y_grid, ymin=ci_l, ymax=ci_u),
                         alpha=1, col=1, width=0.1)
temp_plot = temp_plot + ggplot2::geom_line(data=dataset,
                                           ggplot2::aes(x=y_grid, y=est), lwd=0.8)
temp_plot = temp_plot +
  ggplot2::geom_line(data=dataset, ggplot2::aes(x=y_grid, y=true_dens), col=2, lwd=0.8)
temp_plot = temp_plot +
  ggplot2::geom_ribbon(data=dataset, ggplot2::aes(x=y_grid, ymin=cb_l, ymax=cb_u),
                       alpha=0.2, fill=1)
temp_plot


type = "rbc"
mu = 1
p = 2
true_dens = dnorm(y_grid)

# initializing
est = matrix(0L, nrow=ng, ncol=1)
ci_l = matrix(0L, nrow=ng, ncol=1)
ci_u = cb_l = cb_u = ci_l

# simulate data
x_data = as.matrix(rnorm(5*n))
x_data = x_data[abs(x_data)<=3.2][1:n]
y_data = as.matrix(rnorm(5*n))
y_data = y_data[abs(y_data)<=3.2][1:n]

# estimating bandwidth
bw = lpcde::lpbwcde(y_data=y_data, x_data=x_data, x=x, y_grid = y_grid,
                    p=p, mu=mu, kernel_type="epanechnikov", bw_type="imse-rot")$BW[,2]

# estimating density and confidence bands
model1 = lpcde::lpcde(y_data=y_data, x_data=x_data, x=x, y_grid = y_grid,
                      mu=mu, p=p, bw=bw, kernel_type = "epanechnikov", var_type="ustat")

est = est + model1$Estimate[,3]

# confidence bands
corrMat = sqrt(abs(sweep(sweep(model1$CovMat$CovMat_RBC, MARGIN=1, FUN="*",
                               STATS=1/model1$Estimate[, "se_RBC"]), MARGIN=2, FUN="*",
                         STATS=1/model1$Estimate[, "se_RBC"])))
normalSimu = try(
  MASS::mvrnorm(n=CIsimul, mu=rep(0,nrow(corrMat)),
                Sigma=Matrix::nearPD(corrMat)$mat), silent=TRUE)
z_val = stats::quantile(apply(normalSimu, MARGIN=1,
                              FUN=function(x) {max(abs(x))}), 1 - alpha)
print(z_val)
ci_l = ci_l + (model1$Estimate[,4] - 1.96*sqrt(abs(model1$Estimate[,6])))
ci_u = ci_u + (model1$Estimate[,4] + 1.96*sqrt(abs(model1$Estimate[,6])))
cb_l = cb_l + (model1$Estimate[,4] - z_val*sqrt(abs(model1$Estimate[,6])))
cb_u = cb_u + (model1$Estimate[,4] + z_val*sqrt(abs(model1$Estimate[,6])))

# keeping only alternate confidnece intervals
for (i in 1:ng){
  if (i%%2==0){
    ci_l[i]=NA
    ci_u[i]=NA
  }
}

# creating dataframe
dataset = data.frame(cbind(y_grid, est, true_dens, ci_l, ci_u, cb_l, cb_u))
colnames(dataset) = c("y_grid", "est", "true_dens", "ci_l", "ci_u", "cb_l", "cb_u")

# plotting
temp_plot = ggplot2::ggplot() + ggplot2::theme_bw()
temp_plot = temp_plot+ ggplot2::ylim(-0.5, 0.7)
temp_plot = temp_plot + ggplot2::labs(x="", y="Density") + ggplot2::ggtitle("")
temp_plot = temp_plot +
  ggplot2::geom_errorbar(data=dataset, ggplot2::aes(x=y_grid, ymin=ci_l, ymax=ci_u),
                         alpha=1, col=1, width=0.1)
temp_plot = temp_plot + ggplot2::geom_line(data=dataset,
                                           ggplot2::aes(x=y_grid, y=est), lwd=0.8)
temp_plot = temp_plot +
  ggplot2::geom_line(data=dataset, ggplot2::aes(x=y_grid, y=true_dens), col=2, lwd=0.8)
temp_plot = temp_plot +
  ggplot2::geom_ribbon(data=dataset, ggplot2::aes(x=y_grid, ymin=cb_l, ymax=cb_u),
                       alpha=0.2, fill=1)
temp_plot


#First Derivative Plots
type = "wbc"
mu = 2
p = 3
true_dens = -(y_grid)*exp(-0.5 * (y_grid)^2)/sqrt(2*pi)

# initializing
est = matrix(0L, nrow=ng, ncol=1)
ci_l = matrix(0L, nrow=ng, ncol=1)
ci_u = cb_l = cb_u = ci_l

# simulate data
x_data = as.matrix(rnorm(5*n))
x_data = x_data[abs(x_data)<=3.2][1:n]
y_data = as.matrix(rnorm(5*n))
y_data = y_data[abs(y_data)<=3.2][1:n]

# estimating bandwidth
bw = lpcde::lpbwcde(y_data=y_data, x_data=x_data, x=x, y_grid = y_grid,
                    p=p, mu=mu, kernel_type="epanechnikov", bw_type="imse-rot")$BW[,2]

# estimating density and confidence bands
model1 = lpcde::lpcde(y_data=y_data, x_data=x_data, x=x, y_grid = y_grid,
                      mu=mu, p=p, bw=bw, kernel_type = "epanechnikov", var_type="ustat")

est = est + model1$Estimate[,3]

# confidence bands
corrMat = sqrt(abs(sweep(sweep(model1$CovMat$CovMat, MARGIN=1, FUN="*",
                               STATS=1/model1$Estimate[, "se"]), MARGIN=2, FUN="*",
                         STATS=1/model1$Estimate[, "se"])))
normalSimu = try(
     MASS::mvrnorm(n=CIsimul, mu=rep(0,nrow(corrMat)),
                   Sigma=Matrix::nearPD(corrMat)$mat), silent=TRUE)
z_val = stats::quantile(apply(normalSimu, MARGIN=1,
                              FUN=function(x) {max(abs(x))}), 1 - alpha)
ci_l = ci_l + (model1$Estimate[,3] - 1.96*sqrt(abs(model1$Estimate[,5])))
ci_u = ci_u +(model1$Estimate[,3] + 1.96*sqrt(abs(model1$Estimate[,5])))
cb_l = cb_l + (model1$Estimate[,3] - z_val*sqrt(abs(model1$Estimate[,5])))
cb_u = cb_u + (model1$Estimate[,3] + z_val*sqrt(abs(model1$Estimate[,5])))

# keeping only alternate confidnece intervals
for (i in 1:ng){
  if (i%%2==0){
    ci_l[i]=NA
    ci_u[i]=NA
  }
}

# creating dataframe
dataset = data.frame(cbind(y_grid, est, true_dens, ci_l, ci_u, cb_l, cb_u))
colnames(dataset) = c("y_grid", "est", "true_dens", "ci_l", "ci_u", "cb_l", "cb_u")

# plotting
temp_plot = ggplot2::ggplot() + ggplot2::theme_bw()
temp_plot = temp_plot+ ggplot2::ylim(-1.2, 1.2)
temp_plot = temp_plot + ggplot2::labs(x="", y="Density derivative") + ggplot2::ggtitle("")
temp_plot = temp_plot +
  ggplot2::geom_errorbar(data=dataset, ggplot2::aes(x=y_grid, ymin=ci_l, ymax=ci_u),
                         alpha=1, col=1, width=0.1)
temp_plot = temp_plot + ggplot2::geom_line(data=dataset,
                                           ggplot2::aes(x=y_grid, y=est), lwd=0.8)
temp_plot = temp_plot +
  ggplot2::geom_line(data=dataset, ggplot2::aes(x=y_grid, y=true_dens), col=2, lwd=0.8)
temp_plot = temp_plot +
  ggplot2::geom_ribbon(data=dataset, ggplot2::aes(x=y_grid, ymin=cb_l, ymax=cb_u),
                       alpha=0.2, fill=1)
temp_plot


type = "rbc"
mu = 2
p = 3
true_dens = -(y_grid)*exp(-0.5 * (y_grid)^2)/sqrt(2*pi)

# initializing
est = matrix(0L, nrow=ng, ncol=1)
ci_l = matrix(0L, nrow=ng, ncol=1)
ci_u = cb_l = cb_u = ci_l

# simulate data
x_data = as.matrix(rnorm(5*n))
x_data = x_data[abs(x_data)<=3.2][1:n]
y_data = as.matrix(rnorm(5*n))
y_data = y_data[abs(y_data)<=3.2][1:n]

# estimating bandwidth
bw = lpcde::lpbwcde(y_data=y_data, x_data=x_data, x=x, y_grid = y_grid,
                    p=p, mu=mu, kernel_type="epanechnikov", bw_type="imse-rot")$BW[,2]

# estimating density and confidence bands
model1 = lpcde::lpcde(y_data=y_data, x_data=x_data, x=x, y_grid = y_grid,
                      mu=mu, p=p, bw=bw, kernel_type = "epanechnikov", var_type="ustat")

est = est + model1$Estimate[,3]

# confidence bands
corrMat = sqrt(abs(sweep(sweep(model1$CovMat$CovMat_RBC, MARGIN=1, FUN="*",
                               STATS=1/model1$Estimate[, "se_RBC"]), MARGIN=2, FUN="*",
                         STATS=1/model1$Estimate[, "se_RBC"])))
normalSimu = try(
  MASS::mvrnorm(n=CIsimul, mu=rep(0,nrow(corrMat)),
                Sigma=Matrix::nearPD(corrMat)$mat), silent=TRUE)
z_val = stats::quantile(apply(normalSimu, MARGIN=1,
                              FUN=function(x) {max(abs(x))}), 1 - alpha)
print(z_val)
ci_l = ci_l + (model1$Estimate[,4] - 1.96*sqrt(abs(model1$Estimate[,6])))
ci_u = ci_u + (model1$Estimate[,4] + 1.96*sqrt(abs(model1$Estimate[,6])))
cb_l = cb_l + (model1$Estimate[,4] - z_val*sqrt(abs(model1$Estimate[,6])))
cb_u = cb_u + (model1$Estimate[,4] + z_val*sqrt(abs(model1$Estimate[,6])))

# keeping only alternate confidnece intervals
for (i in 1:ng){
  if (i%%2==0){
    ci_l[i]=NA
    ci_u[i]=NA
  }
}

# creating dataframe
dataset = data.frame(cbind(y_grid, est, true_dens, ci_l, ci_u, cb_l, cb_u))
colnames(dataset) = c("y_grid", "est", "true_dens", "ci_l", "ci_u", "cb_l", "cb_u")

# plotting
temp_plot = ggplot2::ggplot() + ggplot2::theme_bw()
temp_plot = temp_plot+ ggplot2::ylim(-1.2, 1.2)
temp_plot = temp_plot + ggplot2::labs(x="", y="Density derivative") + ggplot2::ggtitle("")
temp_plot = temp_plot +
  ggplot2::geom_errorbar(data=dataset, ggplot2::aes(x=y_grid, ymin=ci_l, ymax=ci_u),
                         alpha=1, col=1, width=0.1)
temp_plot = temp_plot + ggplot2::geom_line(data=dataset,
                                           ggplot2::aes(x=y_grid, y=est), lwd=0.8)
temp_plot = temp_plot +
  ggplot2::geom_line(data=dataset, ggplot2::aes(x=y_grid, y=true_dens), col=2, lwd=0.8)
temp_plot = temp_plot +
  ggplot2::geom_ribbon(data=dataset, ggplot2::aes(x=y_grid, ymin=cb_l, ymax=cb_u),
                       alpha=0.2, fill=1)
temp_plot
