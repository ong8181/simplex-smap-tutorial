####
#### Forecatsint ecological time series using empirical dynamic modeling:
#### a tutorial for simplex projection and S-map
####
#### Functions
####
#### 2018.7.4 Masayuki Ushio & Kazutaka Kawatsu
####

# Define simplex projection and S-map functions
# Implementation of simplex projection and S-map as a function
simplex.projection <- function(ts.d, E){
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
  return(output)
}

s.map <- function(ts.d, E, theta){
  # ts.d = time series data
  # E = number of dimensions for the attractor reconstruction
  # Alternatively, you may use a half of time series to reconstruct attractor and forecast the other half.
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
  return(output)
}


