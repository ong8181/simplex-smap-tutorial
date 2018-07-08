####
#### Forecatsint ecological time series using empirical dynamic modeling:
#### a tutorial for simplex projection and S-map
####
#### 2018.7.8 Masayuki Ushio & Kazutaka Kawatsu


#--------------- Load functions ---------------#
# Load simplex projection and S-map functions (only for Figs.5 & 10)
source("SimplexSmapFuncs.R")

#--------------------------------------------------#
#--------------- Simplex projection ---------------#
#--------------------------------------------------#

#----- Model time series -----#
d <- read.csv("data/TwoSppUni_Example.csv")
tl <- nrow(d)
time <- d$t
# Usually data should be normalized (e.g., using scale() function) before EDM.
# Here we do not perform normalization just for better visualization.
x <- d$X
y <- d$Y


#----- Figure 1 -----#
par(las = 1, mar = c(4,4,3,6), xpd = T)
plot(time, x, type = "l", xlab = "Time", ylab = "Value", lwd = 0.5)
points(time, x, pch = 21, bg = "white", cex = 0.8)
lines(time, y)
points(time, y, pch = 21, bg = gray(0.4), cex = 0.8)
points(time[tl-1], y[tl-1], pch = 23, cex = 1.5, bg = "white")
points(time[tl], y[tl], pch = 23, cex = 1.5, bg = "gray")
legend(par()$usr[2]+1, par()$usr[4], legend = c("X","Y"), lty = 1, pch = c(21, 21), pt.bg = c("white", gray(0.4)))


#----- Embedding -----#
# Assume that the best embedding dimension (E) = 2.
y.E <- 2
y.embed <- matrix(NA, nrow = length(y), ncol = y.E)
for(i in 1:y.E) y.embed[i:length(y), i] <- y[1:(length(y)-i+1)]

# Calculating Euclidean distances between the target and other points
set.target <- matrix(rep(y.embed[tl-1,], nrow(y.embed)), ncol = ncol(y.embed), byrow = T)
distances <- sqrt(rowSums((y.embed - set.target)^2))
neighbors <- order(distances[-(tl-1)])[1:(y.E+1)] # Exclude the target itself


#----- Figure 2 -----#
par(las = 1, mfrow = c(1, 2), mar = c(5,5,2,0.5), xpd = T)
plot(y.embed, xlab = "Y(t)", ylab = "Y(t-1)", pch = 21, bg = "gray", cex = 0.8)
points(y.embed[tl-1,1], y.embed[tl-1,2], pch = 23, bg = "white", cex = 1.5)
for(i in neighbors) points(y.embed[i,1], y.embed[i,2], bg = "black", pch = 23, cex = 1.5)
text(par()$usr[1]-0.03, par()$usr[4]+0.02, "a", cex = 1.3, font = 2)

plot(y.embed, xlab = "Y(t)", ylab = "Y(t-1)", pch = 21, bg = "gray", cex = 0.8)
points(y.embed[tl-1,1], y.embed[tl-1,2], pch = 23, bg = "white", cex = 1.5)
for(i in neighbors){
  arrows(y.embed[i,1], y.embed[i,2],
         y.embed[i+1,1] + 0.008, y.embed[i+1,2] - 0.008, length = 0.1, lwd = 1)
  points(y.embed[i,1], y.embed[i,2], bg = "black", pch = 23, cex = 1.5)
  points(y.embed[i+1,1], y.embed[i+1,2], bg = "gray", pch = 23, cex = 1.5)
}

arrows(y.embed[tl-1,1], y.embed[tl-1,2],
       y.embed[tl,1] + 0.008, y.embed[tl,2] - 0.008, lty = 1, length = 0.1, lwd = 2)
points(y.embed[tl-1,1], y.embed[tl-1,2], pch = 23, bg = "white", cex = 1.5)
points(y.embed[tl,1], y.embed[tl,2], pch = 23, bg = "white", cex = 1.5)
text(par()$usr[1]-0.03, par()$usr[4]+0.02, "b", cex = 1.3, font = 2)


#----- Figure 3 -----#
par(las = 1, mfrow = c(2,1), mar = c(5,5,3,6), xpd = T)
plot(time, y, type = "l", xlab = "Time", ylab = "Value", lwd = 0.5)
points(time, y, pch = 21, bg = "gray", cex = 0.8)
segments(time[tl-1], y[tl-1], d$t[tl-2], d$Y[tl-2], col = "black", lwd = 3)
points(time[tl-1], y[tl-1], pch = 23, cex = 1.5, bg = "white")
points(time[tl-2], y[tl-2], pch = 23, cex = 1.5, bg = "white")
for(i in neighbors){
  points(time[i], y[i], pch = 23, cex = 1.5, bg = "black")
  points(time[i-1], y[i-1], pch = 23, cex = 1.5, bg = "black")
  segments(time[i], y[i], time[i-1], y[i-1], col = "black", lwd = 3)
}
text(par()$usr[1]-5, par()$usr[4]+0.025, "a", cex = 1.3, font = 2)

# Compute weights and calculate prediction
min_distance <- distances[neighbors[1]]
if(min_distance == 0){
  weights <- rep.int(0.000001, times = num_neighbors)
  weights[distances[neighbors] == 0] <- 1
}else{
  weights <- exp(-distances[neighbors]/min_distance)
  weights[weights < 0.000001] <- 0.000001
}
total_weight <- sum(weights)
pred <- (weights %*% y[neighbors+1]) / total_weight

plot(time, y, type = "l", xlab = "Time", ylab = "Value", lwd = 0.5)
points(time, y, pch = 21, bg = "gray", cex = 0.8)

segments(time[tl-1], y[tl-1], time[tl], pred, col = gray(0.4), lwd = 4)
segments(time[tl-1], y[tl-1], d$t[tl-2], d$Y[tl-2], col = "black", lwd = 3)
points(time[tl-1], y[tl-1], pch = 23, cex = 1.5, bg = "white")
points(time[tl-2], y[tl-2], pch = 23, cex = 1.5, bg = "white")
points(time[tl], pred, pch = 23, cex = 1.5, bg = "gray")

for(i in neighbors){
  segments(time[i], y[i], time[i+1], y[i+1], col = gray(0.4), lwd = 4)
  segments(time[i+1], y[i+1], time[tl], y[i+1], lty = 2)
  points(time[i+1], y[i+1], pch = 23, cex = 1.5, bg = "gray")
  points(time[tl], y[i+1], pch = 8, cex = 1.5)
  
  segments(time[i], y[i], time[i-1], y[i-1], col = "black", lwd = 3)
  points(time[i], y[i], pch = 23, cex = 1.5, bg = "black")
  points(time[i-1], y[i-1], pch = 23, cex = 1.5, bg = "black")
}
text(par()$usr[1]-5, par()$usr[4]+0.025, "b", cex = 1.3, font = 2)


#----- Calculating weights -----#
# Compute weights
min_distance <- distances[neighbors[1]]
if(min_distance == 0){
  weights <- rep.int(0.000001, times = num_neighbors)
  weights[distances[neighbors] == 0] <- 1
}else{
  weights <- exp(-distances[neighbors]/min_distance)
  weights[weights < 0.000001] <- 0.000001
}
total_weight <- sum(weights)


#----- Prediction -----#
# Perform one-step prediction
pred <- (weights %*% y[neighbors+1]) / total_weight

# ts.d = time series data
# E = number of dimensions for the attractor reconstruction
# Alternatively, one may use a half of time series to reconstruct attractor.
ts.d <- y
E <- y.E
n <- length(ts.d)

# Make an output dataframe
output <- data.frame(index = 1:n, ts = ts.d, pred = rep(NA, n))

# Univariate SSR
ts.embed <- matrix(NA, nrow = n, ncol = E)
for(i in 1:E) ts.embed[i:n, i] <- ts.d[1:(n-i+1)]

# Calculate distances and make predictions
for(i in 1:(n-1)){
  if(!any(is.na(ts.embed[i,]))){ # skip if the vector contains NA
    set.target <- matrix(rep(ts.embed[i,], n), ncol = E, byrow = T)
    distances <- sqrt(rowSums((ts.embed - set.target)^2))
    neighbors <- order(distances)[-1][1:(E+1)] # Exclude the target itself
    
    # Compute weights
    min.distance <- distances[neighbors[1]]
    if(min.distance == 0){
      weights <- rep.int(0.000001, times = num.neighbors)
      weights[distances[neighbors] == 0] <- 1
    }else{
      weights <- exp(-distances[neighbors]/min.distance)
      weights[weights < 0.000001] <- 0.000001
    }
    total.weight <- sum(weights)
    
    # One-step prediction
    output$pred[i+1] <- (weights %*% ts.d[neighbors+1]) / total.weight
  }
}

# Quantify the forecasting skill by Pearson's correlation coefficient
output <- output[complete.cases(output),]
cor(output$ts, output$pred, method = "pearson")


#----- Figure 4 -----#
par(las = 1, mfrow = c(1,1), xpd = F)
plot(output$ts, output$pred, xlab = "Observed values", ylab = "Predicted values")
abline(0, 1, lty = 2)


#----- Estimating best embedding dimension -----#
# Evaluate forecasting skill for different E
rho.summary <- data.frame(E = 1:10, rho = NA, mae = NA, rmse = NA)

for(E.i in 1:10){
  # Simplex projection is implmented with the function "simplex.projection"
  # Definition can be found Rmarkdown file in GitHub
  simplex.res <- simplex.projection(y, E.i)
  simplex.res <- simplex.res[complete.cases(simplex.res),]
  rho.summary$rho[E.i] <- cor(simplex.res$ts, simplex.res$pred, method = "pearson")
  rho.summary$mae[E.i] <- mean(abs(simplex.res$pred - simplex.res$ts))
  rho.summary$rmse[E.i] <- sqrt(mean((simplex.res$pred - simplex.res$ts)^2))
}

#----- Figure 5 -----#
par(las = 1, mfrow = c(1,3), mar = c(5,5,2,0.5), xpd = T)
plot(rho.summary$E, rho.summary$rho, type = "b", xlab = "E", ylab = expression(rho))
text(par()$usr[1]-0.8, par()$usr[4]+0.02, "a", cex = 1.5, font = 2)
plot(rho.summary$E, rho.summary$mae, type = "b", xlab = "E", ylab = "MAE")
text(par()$usr[1]-0.8, par()$usr[4]+0.004, "b", cex = 1.5, font = 2)
plot(rho.summary$E, rho.summary$rmse, type = "b", xlab = "E", ylab = "RMSE")
text(par()$usr[1]-0.8, par()$usr[4]+0.004, "c", cex = 1.5, font = 2)


#-------------------------------------------------#
#--------------------- S-map ---------------------#
#-------------------------------------------------#

#----- Fitting a locally weighted linear model -----#
#----- Figure 6 -----#
x.grid <- seq(min(y.embed[,1], na.rm=T)-0.05, max(y.embed[,1], na.rm=T)+0.05, by=0.015)
y.grid <- seq(min(y.embed[,2], na.rm=T)-0.05, max(y.embed[,2], na.rm=T)+0.05, by=0.015)
x.axis <- rep(x.grid, length(y.grid))
y.axis <- NULL
for(i in 1:length(y.grid)) y.axis <- c(y.axis, rep(y.grid[i], length(x.grid)))
grid.xy <- data.frame(cbind(x.axis,y.axis))
set.target2 <- matrix(rep(y.embed[tl-1,], nrow(grid.xy)), ncol = ncol(y.embed), byrow = T)

grid.xy$dist <- sqrt(rowSums((grid.xy - set.target2)^2))
grid.ws <- exp(-2 * grid.xy$dist/mean(distances, na.rm = T)) # Theta is set as 2 in this case
grid.ws.mat <- matrix(grid.ws, ncol = length(x.grid))

par(las = 1, mfrow = c(1,1), xpd = F)
plot(y.embed[,1], y.embed[,2], pch = 21, cex = 0, xlab = "Y(t)", ylab = "Y(t-1)")
image(x.grid, y.grid, grid.ws.mat, col = gray.colors(100, start = 1, end = 0.3), add = T)
points(y.embed[,1], y.embed[,2], pch = 21, bg = "gray", cex = 0.8)
points(y.embed[tl-1,1], y.embed[tl-1,2], pch = 23, cex = 1.5, bg = "white")


#----- Figure 7 -----#
# Weights for all points
plot.ws <- exp(-2 * distances/mean(distances, na.rm = T)) # Theta is set as 2 in this case

par(las = 1, mar = c(4,4,3,6), xpd = T)
plot(time, y, type = "l", lwd = 0.5, xlab = "Time", ylab = "Value")
# Add segment color
for(i in 2:length(y)){
  if(i != tl-1) segments(time[i], y[i], time[i-1], y[i-1], col = gray(1 - plot.ws[i]), lwd = 5*plot.ws[i])
}
points(time, y, pch = 21, bg = "gray")

segments(time[tl-1], y[tl-1], d$t[tl-2], d$Y[tl-2], col = "black")
points(time[tl-1], y[tl-1], pch = 23, cex = 1.5, bg = "white")
points(time[tl-2], y[tl-2], pch = 23, cex = 1.5, bg = "white")


#----- Figure 8 -----#
# ts.d = time series data
# E = number of dimensions for the attractor reconstruction
# Alternatively, one may use a half of time series to reconstruct attractor.
ts.d <- y
E <- y.E
n <- length(ts.d)

# Univariate SSR
ts.embed <- matrix(NA, nrow = n, ncol = E)
for(i in 1:E) ts.embed[i:n, i] <- ts.d[1:(n-i+1)]

# Calculate distances and make predictions
set.target <- matrix(rep(ts.embed[tl-1,], n), ncol = E, byrow = T)
distances <- sqrt(rowSums((ts.embed - set.target)^2))

# Compute weights
d.m <- mean(distances, na.rm = T)
ws <- exp(-10 * distances / d.m) # Set theta = 10 here

# Make data.frame w/o NA and target vector
ts.ws0 <- cbind(ts.embed, ws, matrix(rep(NA, n), ncol = 1))
ts.ws0[1:(n-1),4] <- ts.embed[2:n,1]
ts.ws <- ts.ws0[-(tl-1),]
ts.ws <- ts.ws[complete.cases(ts.ws),]
ts.ws1 <- ts.ws
ts.ws1[,3] <- 1 # For unweighted linear regression

# Weighted linear regression (S-map)
# Singular-value decomposition (excluding the target vector itself)
A <- cbind(ts.ws[,1:2], 1) * ts.ws[,3]
A.svd <- svd(A)
# Remove singular values that are too small
s <- A.svd$d
s.inv <- matrix(0, nrow = E+1, ncol = E+1)
for(j in seq_along(s)) s.inv[j,j] <- 1/s[j]
# Perform back-substitute to solve        
map <- A.svd$v %*% s.inv %*% t(A.svd$u) %*% (ts.ws[,3] * ts.ws[,4])

# Unweighted linear regression
# Singular-value decomposition (excluding the target vector itself)
A1 <- cbind(ts.ws1[,1:2], 1) * ts.ws1[,3]
A1.svd <- svd(A1)
# Remove singular values that are too small
s1 <- A1.svd$d
s1.inv <- matrix(0, nrow = E+1, ncol = E+1)
for(j in seq_along(s1)) s1.inv[j,j] <- 1/s1[j]
# Perform back-substitute to solve        
map1 <- A1.svd$v %*% s1.inv %*% t(A1.svd$u) %*% (ts.ws1[,3] * ts.ws1[,4])

# Make prediction
output0 <- cbind(ts.embed, 1) %*% map
output1 <- cbind(ts.embed, 1) %*% map1

smap.resid0 <- (output0[1:169] - ts.d[2:170])
smap.resid1 <- (output1[1:169] - ts.d[2:170])
out.df <- data.frame(dist = distances[1:169], weighte_resid = smap.resid0, unweighted_resid =smap.resid1)

# Visualization
par(las = 1, mar = c(4,4,3,6), xpd = T)
plot(out.df[,1], out.df[,2], pch = 21, bg = gray(0.5),
     ylim = c(-0.8, 0.2), xlab = "Distance from the target vector", ylab = "S-map residual")
points(out.df[,1], out.df[,3], pch = 21, bg = "white")
segments(-0.01, 0, 0.365, 0, lty = 2)
legend(par()$usr[2]+0.01, par()$usr[4],
       legend = c(expression(paste(theta, " = 0")),
                  expression(paste(theta, " = 10"))),
       pch = c(21, 21), pt.bg = c("white", gray(0.5)))


#----- Implementation of S-map with R -----#
# ts.d = time series data
# E = number of dimensions for the attractor reconstruction
# Alternatively, one may use a half of time series to reconstruct attractor.
ts.d <- y
E <- y.E
n <- length(ts.d)
theta <- 10

# Make an output dataframe
output <- data.frame(index = 1:n, ts = ts.d, pred = rep(NA, n))

# Univariate SSR
ts.embed <- matrix(NA, nrow = n, ncol = E)
for(i in 1:E) ts.embed[i:n, i] <- ts.d[1:(n-i+1)]

# Calculate distances and make predictions
for(i in 1:(n-1)){
  if(!any(is.na(ts.embed[i,]))){ # skip if the vector contains NA
    set.target <- matrix(rep(ts.embed[i,], n), ncol = E, byrow = T)
    distances <- sqrt(rowSums((ts.embed - set.target)^2))
    
    # Compute weights
    d.m <- mean(distances, na.rm = T)
    ws <- exp(-theta * distances / d.m)
    
    # Make data frame w/o NA and target vector
    ts.ws0 <- cbind(ts.embed, ws, matrix(rep(NA, n), ncol = 1))
    ts.ws0[1:(n-1),4] <- ts.embed[2:n,1]
    ts.ws <- ts.ws0[-i,]
    ts.ws <- ts.ws[complete.cases(ts.ws),]
    
    # Singular-value decomposition (excluding the target vector itself)
    A <- cbind(ts.ws[,1:2], 1) * ts.ws[,3]
    A.svd <- svd(A)
    
    # Remove singular values that are too small
    s <- A.svd$d
    s.inv <- matrix(0, nrow = E+1, ncol = E+1)
    for(j in seq_along(s))
    {
      if(s[j] >= max(s) * 1e-5)
        s.inv[j,j] <- 1/s[j]
    }
    
    # Perform back-substitute to solve        
    map <- A.svd$v %*% s.inv %*% t(A.svd$u) %*% (ts.ws[,3] * ts.ws[,4])
    
    # Make prediction
    output$pred[i+1] <- sum(map * c(ts.embed[i,], 1))
  }
}

# Quantify the forecasting skill
output <- output[complete.cases(output),]
cor(output$ts, output$pred, method = "pearson")


#----- Figure 9 -----#
par(las = 1, mfrow = c(1,1), xpd = F)
plot(output$ts, output$pred, xlab = "Observed values", ylab = "Predicted values")
abline(0, 1, lty = 2)


#----- Quantifying nonlinearlity -----#
# Evaluate forecasting skill for different theta
smap.summary <- data.frame(theta = 0:10, rho = NA, mae = NA, rmse = NA)

for(theta.i in 0:10){
  # S-map is implmented with the function "s.map"
  # Definition can be found Rmarkdown file in GitHub
  smap.res <- s.map(y, 2, theta.i)
  smap.res <- smap.res[complete.cases(smap.res),]
  smap.summary$rho[theta.i+1] <- cor(smap.res$ts, smap.res$pred, method = "pearson")
  smap.summary$mae[theta.i+1] <- mean(abs(smap.res$pred - smap.res$ts))
  smap.summary$rmse[theta.i+1] <- sqrt(mean((smap.res$pred - smap.res$ts)^2))
}


#----- Figure 10 -----#
par(las = 1, mfrow = c(1,3), mar = c(5,5,2,0.5), xpd = T)
plot(smap.summary$theta, smap.summary$rho, type = "b", xlab = expression(theta), ylab = expression(rho))
text(par()$usr[1]-0.8, par()$usr[4]+0.02, "a", cex = 1.5, font = 2)
plot(smap.summary$theta, smap.summary$mae, type = "b", xlab = expression(theta), ylab = "MAE")
text(par()$usr[1]-0.8, par()$usr[4]+0.004, "b", cex = 1.5, font = 2)
plot(smap.summary$theta, smap.summary$rmse, type = "b", xlab = expression(theta), ylab = "RMSE")
text(par()$usr[1]-0.8, par()$usr[4]+0.004, "c", cex = 1.5, font = 2)

