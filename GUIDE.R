library(dplyr)
library(TH.data)
library(prodlim)
library(MASS)

# http://publicifsv.sund.ku.dk/~tag/Teaching/share/R-tutorials/Advanced-statistics/SurvivalAnalysis.html
# https://www.rdocumentation.org/packages/TH.data/versions/1.0-10/topics/GBSG2
data("GBSG2")

# X - K dimensions vector of covariates
# Y - univariate response variable 
# Z - treatment variable 1 ... G

Y <- GBSG2[[ncol(GBSG2)]]
Z <- GBSG2[[1]]
X <- GBSG2[2:(ncol(GBSG2)-1)]

G <- Z %>% levels() %>% length()
K <- X %>% ncol()

# add some NA values
# TODO that functions will work with nans 
if(F){
  n <- as.integer( nrow(X)*0.05 )
  c <- sample(1:ncol(X), n, replace=T)
  r <- sample(1:nrow(X), n, replace=T)
  for(i in 1:n){
    X[r[i],c[i]] <- NA
  }
  # add NA as factor
  for(i in 1:K){
    if( is.ordered(X[[i]])){
      X[[i]] <- factor(X[[i]], levels = c(NA, levels(X[[i]])), exclude = "", ordered = T)
    }else if(is.factor(X[[i]])){
      X[[i]] <- addNA(X[[i]],ifany = FALSE)
    }
  }
}

summary(X)
summary(Y)
summary(Z)

split_variable_selection <- function(Xj,Y,Z,G){
  # input: X[[]] (one column)
  # 1
  if(is.factor(Xj)){
    V <- Xj %>% as.numeric()
    h <- levels(Xj) %>% length()
  # 2
  }else{
    # a
    m <- Xj %>% unique() %>% length()
    has.na <- Xj %>% is.na() %>% sum() > 0
    nt <- Xj %>% length()

    if(m<5||(m==5&&has.na)){ 
      h <- m
    }else if (nt < 30*(G+1)){ 
      h <- 3
    }else{ 
      h <- 4
    }
    
    if(has.na){ 
      rk <- sapply(1:(h-2), function(k) k/(h-1))
    }else{ 
      rk <- sapply(1:(h-1), function(k) k/h) 
    }
    #b 
    q0 <- -Inf
    qk <- sapply(rk, function(r) Xj[floor(r*nt)])
    qk <- c(q0, qk)
    
    if(has.na){ 
      V <- sapply(1:(h-2), 
                  function(k) k*(Xj>qk[k] & Xj<=qk[k+1] & !is.na(Xj))) %>% rowSums()
      V <- V + (h-1)*(Xj>qk[h-2+1] & !is.na(Xj)) + h*is.na(Xj)
    }else{ 
      V <- sapply(1:(h-1), 
                  function(k) k*(Xj>qk[k]&Xj<=qk[k+1])) %>% rowSums()
      V <- V + (h-1)*(Xj>qk[h-1+1]) }
  }
  # 3
  # TEST
  additive_data <- data.frame(y=Y)
  for(z in 1:G){
    additive_data[paste0('z',z)] <- as.numeric(Z==z) 
  }
  for(v in 1:length(unique(V))){
    additive_data[paste0('v',v)] <- as.numeric(V==v) 
  }
                             
  add.lm <- lm(y ~ ., data=additive_data)
  
  full_data <- data.frame(y=Y)
  for(z in 1:G){
    for(v in 1:length(unique(V))){
      full_data[paste0('z',z,'v',v)] <- as.numeric((Z==z)&(V==v)) 
    }
  }
  
  full.lm <- lm(y ~ ., data=full_data)
    
  pj <- anova(add.lm, full.lm)$`Pr(>F)`[2]
  return(pj)
}

left_right_c <- function(Xj, Xj_points, j, split.type=2){
  c <- (Xj_points[j+1]+Xj_points[j])/2

  if(split.type==1){
    left <- (is.na(Xj)) %>% which()
    right <- (!is.na(Xj)) %>% which()
  }else if(split.type==2){
    left <- (Xj<=c | is.na(Xj)) %>% which()
    right <- (Xj>c & !is.na(Xj)) %>% which()
  }else if(split.type==3){
    left <- (Xj<=c & !is.na(Xj)) %>% which()
    right <- (Xj>c | is.na(Xj)) %>% which()
  }
  return(list("c"=c, "left"=left, "right"=right))
}

split_set_selection <- function(Xj,Y,Z,G,split.type=2){
  # input: X[[]] (one column)
  if(is.numeric(Xj)){
    # c <- midpoint between consecutive order statistics of X in t (mid points between ordered X values)
    # X has m order statistics -> max number of splits is (m-1) or {1+2(m-1)} <- if NA or not
    #  Permissible splits are those that yield two child nodes with each having two or more
    #  observations per treatment. The selected split is the one that minimizes the sum of
    #  the deviances (or sum of squared residuals in the case of least-squares regression)
    #  in the two child nodes.
    min_sum <- Inf 
    min_j <- 0
    
    Xj_points <- Xj %>% unique() %>% sort()
    # iterate over all possible midpoints
    for(j in 2:(length(Xj_points)-2)){

      list_lrc <- left_right_c(Xj,Xj_points, j, split.type)
      c <- list_lrc$c
      left <- list_lrc$left
      right <- list_lrc$right
      
      # check if min 2 values per treatment
      left_min <- Y[left] %>% table() %>% min()
      right_min <- Y[right] %>% table() %>% min()
      if((left_min > 1)&(left_min < Inf)&(right_min > 1)&(right_min < Inf)){
        # linear model for left nd right
        fit_l <- lm(Y[left] ~ Xj[left])#, na.action=na.omit)
        fit_r <- lm(Y[right] ~ Xj[right])#, na.action=na.omit)
        # coefficients
        l_coeff <- fit_l$coefficients
        r_coeff <- fit_r$coefficients
        # sum of least squares
        l_sum <- sum((Y[left]  - (l_coeff[1] + l_coeff[2]*Xj[left])) ^ 2)
        r_sum <- sum((Y[right]  - (r_coeff[1] + r_coeff[2]*Xj[right])) ^ 2)
        # check sum of least squares if new are better than old
        if(l_sum+r_sum < min_sum){
          min_j <- j
          min_sum <-l_sum+r_sum
        }
      }
    }
    print(min_j)
    list_lrc <- left_right_c(Xj, Xj_points, min_j, split.type)
    left <- list_lrc$left
    right <- list_lrc$right
    # 2
  }else{
    m <- levels(Xj) %>% length()
    if(F && m<=11){
      # for small data testing the second way
      # complete search
      levels_m <- levels(Xj) %>% as.vector()
      print(levels_m)
    
    }else{
      # performs an approximate search by means of linear discriminant analysis
      # 1
      y_z <- lapply(1:G, function(z) Y[as.numeric(Z)==z] %>% mean() )
      # 2
      C <- sapply(1:length(Y), function(y) if(Y[y]>y_z[as.numeric(Z[y])]) 2*as.numeric(Z[y])-1 else 2*as.numeric(Z[y]) )
      # 3
      D <- sapply(1:m, function(a) (as.numeric(Xj)==a) %>% as.numeric() ) %>% data.frame()
      # 4 
      DC <- D
      DC['class'] <- C
      lda_dc <- lda(class ~ ., DC)
      # I am choosing first LDA variable
      Bvals <- lda_dc$scaling[,1]
      B <- (D * lda_dc$scaling[,1]) %>% rowSums()
      # find best way of splitting
      min_sum <- Inf 
      min_b <- 0
      for(b in Bvals[1:(length(Bvals)-1)]){
        left <- (B<=b) %>% which()
        right <- (B>b) %>% which()
        # check if min 2 values per treatment
        left_min <- Y[left] %>% table() %>% min()
        right_min <- Y[right] %>% table() %>% min()
        if((left_min > 1)&(left_min < Inf)&(right_min > 1)&(right_min < Inf)){
          # linear model for left nd right
          print(Y[left] %>% table() )
          print(Y[right] %>% table() %>% min()  )
          fit_l <- lm(Y[left] ~ Xj[left])
          fit_r <- lm(Y[right] ~ Xj[right])
          print(8)
          # coefficients
          l_coeff <- fit_l$coefficients
          r_coeff <- fit_r$coefficients
          # sum of least squares
          l_sum <- sum((Y[left]  - (l_coeff[1] + l_coeff[2]*Xj[left])) ^ 2)
          r_sum <- sum((Y[right]  - (r_coeff[1] + r_coeff[2]*Xj[right])) ^ 2)
          print(l_coeff)
          print(r_sum)
          # check sum of least squares if new are better than old
          if((l_sum+r_sum) < min_sum){
            min_b <- b
            min_sum <-l_sum+r_sum
          }
        }
      }
      left <- (B<=min_b) %>% which()
      right <- (B>min_b) %>% which()
    }
  }
  return(list("left"=left, "right"=right))
}