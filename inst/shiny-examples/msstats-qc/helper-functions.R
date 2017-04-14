getPolarCoord <- function(r, matrix = F, na = F){
  # Get starting angle and angle increments
  theta <- 0
  dtheta <- 360 / length(r)
  dtheta <- (pi / 180) * dtheta  # in radians
  
  # Get polar coordinates
  x <- c()
  y <- c()
  
  for(i in 1:length(r)){
    
    x <- c(x, r[i] * cos(theta))
    y <- c(y, r[i] * sin(theta))
    
    theta <- theta + dtheta
  }
  
  x[length(x) + 1] <- x[1]
  y[length(y) + 1] <- y[1]
  
  if(na == T){
    x[length(x) + 1] <- NA
    y[length(y) + 1] <- NA
  }
  
  
  if(matrix == T){
    return(cbind(x, y))
  }else{
    return(list(x = x, 
                y = y))
  }
  
}