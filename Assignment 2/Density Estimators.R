# The implementation of the histogram estimator
Histogram = function(x, h = 1, data){
  max = max(data)
  min = min(data)
  # The default value of x_0 is the mean of the data
  x0 = mean(data)
  breaks = unique(c(tail(seq(x0, min, -h), 1) - h, seq(x0, min, -h),
                    seq(x0, max, h), tail(seq(x0, max, h), 1) + h))
  density = hist(data, breaks = breaks, right = FALSE, plot = FALSE)$density
  n = length(breaks)
  # Give the ability of tackling vectorized input to this function
  p = length(x)
  return = rep(0, p)
  for(i in 1:p){
    location = sum(x[i] >= breaks)
    if((location > 0) & (x[i] <= tail(breaks, 1))){
      return[i] = density[location]
    }
  }
  return
}

# The implementation of the orthogonal series estimator, using Hermite series
OS_Hermite = function(x, K = 30, data){
  # The Hermite polynomials
  # Utilize the recurrence formula "H_n(x) = 2x*H_(n-1)(x) - 2(n-1)*H_(n-2)(x) (n>=2)"
  # to calculate the Hermite polynomials
  H_m = function(x, m){
    if(m == 0){
      1
    }
    else if(m == 1){
      2*x
    }
    else{
      a = 1
      b = 2*x
      for(k in 2:m){
        c = 2*x*b - 2*(k - 1)*a
        a = b
        b = c
      }
      c
    }
  }
  # The normalized Hermite functions
  phi_m = function(x, m){
    H_m(x, m)/exp((x^2)/2)/sqrt((2^m)*factorial(m)*sqrt(pi))
  }
  n = length(data)
  value = 0
  for(m in 0:K){
    c_m = sum(sapply(data, phi_m, m = m))/n
    value = value + c_m*phi_m(x, m)
  }
  value
}

# The implementation of the orthogonal series estimator, using trigonometric series
OS_trigo = function(x, K = 30, data){
  phi_r = function(x, r){
    if(r == 0){
      1
    }
    else if(r%%2 == 0){
      sqrt(2)*sin(pi*r*x)
    }
    else if(r%%2 == 1){
      sqrt(2)*cos(pi*(r + 1)*x)
    }
  }
  n = length(data)
  value = 0
  for(r in 0:K){
    c_r = sum(sapply(data, phi_r, r = r))/n
    value = value + c_r*phi_r(x, r)
  }
  value
}
