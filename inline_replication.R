# Replication file for inline code snippets in Cattaneo, Chandak, Janson and Ma (2022)
# Submitted to the R Journal

#Install package
install.packages("lpcde")
# load package
library("lpcde")

# software article script
set.seed(42)
n=1000
x_data = as.matrix(rnorm(n, mean=0, sd=1))
y_data = as.matrix(rnorm(n, mean=0, sd=1))
y_grid = stats::quantile(y_data, seq(from=0.1, to=0.9, by=0.1))

# density estimation
model1 = lpcde::lpcde(x_data=x_data, y_data=y_data, y_grid=y_grid, x=0, bw=0.5)
summary(model1)

# simple plot
model1 = lpcde::lpcde(x_data=x_data, y_data=y_data, x=0, bw_type = "mse-rot")
plot(model1) + ggplot2::theme(legend.position="none")

# bandiwdth selection
y_grid = stats::quantile(y_data, seq(from=0.1, to=0.9, by=0.1))
model2 = lpcde::lpbwcde(y_data=y_data, x_data=x_data, x=0, y_grid = y_grid, bw_type = "mse-rot")
summary(model2)
