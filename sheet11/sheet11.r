# some part of the code is based on the ideas and/or the implementation of the qsimularR vignettes (https://cran.r-project.org/web/packages/qsimulatR/vignettes/)
library("qsimulatR")


generate_basis <- function(n) {
  # generates a basis for the quantum state with the necessary registers for the 
  # order finding algorithm
  basis <- c()
  total_states <- 2^(2 * n + 4 + (2 * n + 3))   
  for (i in 0:(total_states - 1)) {
    t <- i %/% (2^(2 * n + 4))  
    xbits <- (i %/% (2^(n + 4))) %% 2^n  
    reg2 <- (i %/% (2^4)) %% 2^n
    helpers <- c((i %/% 16) %% 2, (i %/% 8) %% 2, (i %/% 4) %% 2, (i %/% 2) %% 2) 
    
    basis[i + 1] <- paste0(
      "|t=", t, 
      ">|xbits=", xbits, 
      ">|reg2=", reg2, 
      ">|helpers=", paste(helpers, collapse = ""), ">"
    )
  }
  
  return(basis)
}

write_measurement <- function(object, file, ...) {
  # helper function for writing the measurement results to a file due to 
  # to long calculation times

  # Open a connection to the file
  con <- file(file, "w")
  
  if (is.na(object$bit)) {
    cat("All bits have been measured", object$repetitions, "times with the outcome:\n", file = con)
    tmp.state <- qstate(object$nbits, basis = object$basis)
    tmp.state@coefs <- as.complex(object$value)
    capture.output(show(tmp.state), file = con) 
  } else {
    cat("Bit", object$bit, "has been measured", object$repetitions, "times with the outcome:\n", file = con)
    ones <- sum(object$value)
    zeros <- object$repetitions - ones
    cat("0: ", zeros, "\n1: ", ones, "\n", file = con)
  }

  close(con)
}

import_t_values <- function(file_path, threshold){
  # reads the exported values from the function "write_measurement" and extracts the t values
  # into an array
  # threshold is the necessary number of measurements that a value is considered valid and not as random
  lines <- readLines(file_path)

  t_values <- integer()
  
  for (i in 2:length(lines)) {
    line <- lines[i]
    num_measured <- as.integer(sub("^.*\\(\\s*(\\d+)\\s*\\)\\s.*", "\\1", line)) 
    if (!is.na(num_measured) && num_measured > threshold) {
    t_value <- as.integer(sub(".*\\|t=([0-9]+)>.*", "\\1", line))
    t_values <- c(t_values, t_value)
    } 
  }
   
  return(unique(t_values))
}

cadd <- function(c, psi, y, xbits){
  # Add the integer y to the register xbits (representing x) (mod 2^n) where n is the number of bits in x
  # psi has to be a qstate object with c is an additional control bit

  n <- length(xbits)
  q <- cqft(c, psi, bits=xbits)
  for (j in 1:n) {
    theta <- 2 * pi * y / 2^(1-j+n)
    phasegate <- sqgate(bit=j, M=array(as.complex(c(1,0,0,exp(1i*theta) )), dim=c(2,2))) 
    q <- cqgate(bits=c(c,xbits[j]), gate=phasegate) * q
  }
  q <- cqft(c,q, bits=xbits, inverse=TRUE)
  return(q)
}

summands <- function(y, N, n){
  yb <- as.integer(intToBits(y))
  r <- c()
  for (i in c(1:n)) {
    s <- 0
    for (j in c(1:n)) {
      s <- s + as.integer(yb[j]) * 2^(i+j-2)
    }
    r[i] <- s %% N
  }
  return(r)
}

eEa <- function(a, b){
  # implements the extended euclidean algorithm
  if (a == 0) {
    return(c(b, 0, 1))
  }
  res <- eEa(b %% a, a)
  return(c(res[1], res[3] - (b %/% a) * res[2], res[2]))
}

modinv <- function(a, n){
  # calculate the modular inverse of a mod n
  res <- eEa(a, n)
  if (res[1] != 1) {
    stop("a and 2^n are not coprime, the modular inverse does not exist")
  }
  return(res[2] %% n)
}

add <- function(psi, y, xbits){
  # Add the integer y to the register xbits (representing x) (mod 2^n) where n is the number of bits in x
  # psi has to be a qstate object with c is an additional control bit

  n <- length(xbits)
  q <- qft(psi, bits=xbits)
  for (j in 1:n) {
    theta <- 2 * pi * y / 2^(1-j+n)
    phasegate <- sqgate(bit=xbits[j], M=array(as.complex(c(1,0,0,exp(1i*theta) )), dim=c(2,2))) 
    q <- phasegate * q
  }
  q <- qft(q, bits=xbits, inverse=TRUE)
  return(q)
}

is.less <- function(psi, y, a, c1, xbits){
  
  bits <- c(xbits, a) # add the helper bit for the underflow
  n <- length(bits)

  q <- add( psi, y=2^n - y, xbits=bits) # subtract y 
  q <- CNOT(c(a, c1)) * q # store the value of the underflow to the helper qubit c1 

  q <- add(q, y, bits) # add y back to restore the original state beside the value of c1

  return(q)
}

cis.less <- function(c, psi, y, a, c1, xbits){
  
  bits <- c(xbits, a) # add the helper bit for the underflow
  n <- length(bits)

  q <- cadd(c, psi, y=2^n - y, xbits=bits) # subtract y 
  q <- CNOT(c(a, c1)) * q # store the value of the underflow to the helper qubit c1 

  q <- cadd(c, q, y, bits) # add y back to restore the original state beside the value of c1

  return(q)
}

caddmodN <- function(c, psi, y, N, a, c1, c2, xbits){
  ###############
  ###  TASK 1 ###
  ###############
  n <- length(xbits)
  
  # make sure that y is smaller than N
  y <- y %% N

  q <- cis.less(c, psi, y=N, a=a, c1=c1, xbits=xbits) #check if x < N and store the result in c1
  q <- cis.less(c, q, y=N-y, a=a, c1=c2, xbits=xbits) #check if x < N-y and store the result in c2


    
  # if c1=1 and c2=0 compute x=x+y-N
  q <- X(c2) * (CCNOT(c(c1, c2, a))* (X(c2)*q))
  q <- cadd(a, q, y=y-N, xbits=xbits) 
  q <- X(c2)* (CCNOT(c(c1, c2, a))*(X(c2)*q))


  # if c1=1 and c2=1 compute x=x+y
  q <- CCNOT(c(c1, c2, a))*q
  q <- cadd(a, q, y=y, xbits=xbits)
  q <- CCNOT(c(c1, c2, a))*q
  
  #reset helper bits
  q <- cis.less(c, q, y=y, a=a, c1=c2, xbits=xbits)
  q <- CNOT(c(c1, c2)) * q
  q <- cis.less(c, q, y=N, a=a, c1=c1, xbits=xbits)
  
  return(q)
}


cmultmodN <- function(c, psi, y, N, xbits, reg2, helpers){
  ################# 
  ###  TASK 2  ####
  #################
  # controlled version of |x>|0> -> |x*y mod N>|0>
  # xbits has to have the same length as reg1
  # 4 helper bits are required
  n <- length(xbits) # get length of |x>  
  n2 <- length(reg2) # get length of |0>
  stopifnot(n == n2) # check if |x> and |0> have the same length
  s <- summands(y, N, n) # precompute the summands from y
  # calculate x*y mod N and store it in reg2
  for (i in c(1:n)) {
    psi <- CCNOT(c(c, xbits[i], helpers[4])) * psi
    psi <- caddmodN(c=helpers[4], psi, y=s[i], N=N, a=helpers[3], c1=helpers[1], c2=helpers[2], xbits=reg2)
    psi <- CCNOT(c(c, xbits[i], helpers[4])) * psi
  }
  # swap the registers
  for (i in c(1:n)) {
    psi <- CSWAP(c(c, xbits[i], reg2[i])) * psi
  }
  # reset reg2
  yinv <- N-modinv(y, N)
  s <- summands(yinv, N, n)

  for (i in c(1:n)) {
    psi <- CCNOT(c(c, xbits[i], helpers[4])) * psi
    psi <- caddmodN(c=helpers[4], psi, y=s[i], N=N, a=helpers[3], c1=helpers[1], c2=helpers[2], xbits=reg2)
    psi <- CCNOT(c(c, xbits[i], helpers[4])) * psi
  }
  return(psi)
}

cexpmodNs <- function(c, a, psi, y, N, xbits, reg2, helpers){
  # controlled version of |x>|0> -> |x^a*y mod N>|0>
  # xbits has to have the same length as reg1
  # 4 helper bits are required
  for (i in c(1:a)) {
    psi <- cmultmodN(c, psi, y, N, xbits, reg2, helpers)
  }
  return(psi)
}

cexpmodN <- function(c, a, psi, y, N, xbits, reg2, helpers){
  # controlled version of |x>|0> -> |x^a mod N>|0>
  # xbits has to have the same length as reg1
  # 4 helper bits are required
  abin <- as.integer(intToBits(a))
  n <- max(which(abin == 1))
  y2 <- y %% N
for (i in c(1:n)) {
    if (abin[i] == 1) { # if statement to save one bit like qsimulatR vignettes since the code is already very slow in R. If that is not okay, the cexpmodNs function can be used instead
    psi <- cmultmodN(c, psi, y=y2, N=N, xbits=xbits, reg2=reg2, helpers=helpers)
    }
    y2 <- (y2 * y2) %% N
  }
  return(psi)
}

continued_fraction <- function(x, eps=1e-14, k_max=100){
  #################
  ###  TASK 3  ####
  #################
  # compute the continued fraction of x up to a precision of eps
  # or a maximum of k_max iterations
  k <- c()
  p <- c(0, 1)
  q <- c(1, 0)
  s <- c()

  x0 <- x
  for (i in 1:k_max) {
    kn <- floor(x)
    pn <- kn * rev(p)[1] + rev(p)[2]
    qn <- kn * rev(q)[1] + rev(q)[2]
    sn <- pn/qn

    p <- c(p, pn)
    q <- c(q, qn)
    k <- c(k, kn)
    s <- c(s, sn)

    frac_part <- x - kn

     # Break condition if precision is reached
      if (abs(x0 - sn) < eps) {
        break
      }

    x <- 1/frac_part
}
  return(q)
}

qpe_order_finding <- function(t_bits, psi, y, N, xbits, reg2, helpers){
  # apply the phase estimation algorithm to the order finding problem
  psi <- X(xbits[1]) * psi # set the first bit of x to 1 such that is an eigenstate of U 

  for (i in t_bits){
    psi <- H(i) * psi # apply Hadamard gate to the t register
  }

  j <- 1
  for (i in t_bits){ # apply the controlled U gates
    cat(sprintf("Progress: qpe bit %d of %d\n", j, length(t_bits)))
    psi <- cexpmodN(c=i, a=2^(j-1), psi=psi, y=y, N=N, xbits=xbits, reg2=reg2, helpers=helpers)
    j <- j + 1
  }
  psi <- qft(psi, bits=t_bits, inverse=TRUE)
  return(psi)
}

###############################  MAIN  ####################################

# define basis and parameters
n <- 3
basis <- generate_basis(n)
bits_total <- 4*n+7  
tbits_total <- 2*n+3

q <- qstate(bits_total, basis=basis)

helpers <- c(1:4)
reg2 <- c(5:7)
xbits <- c(8:10)
N <- 7
y <- 3
t_bits <- c(11:bits_total)

##############################################################################
# for saving the time of the phase estimation this part can be commented out 
q <- qpe_order_finding(t_bits=t_bits, psi=q, y=y, N=N, xbits=xbits, reg2=reg2, helpers=helpers) # perform the phase estimation algorithm
res <- measure(q, repetitions = 10000) # measure the qubit state
write_measurement(res, "measurement_results.txt") # write values to file
##############################################################################

t_values <- import_t_values("measurement_results.txt", threshold=100) # read all t values which were more than 100 times measured
  
fracts <- c()
for (value in t_values){
  if (value != 0){ # skip the 0 values since it is no useful value for the continued fraction
  frac = continued_fraction(value/2^tbits_total) # convert the t values to float numbers and apply the continued fraction algorithm to get the q from phi= p/q
  fracts <- c(fracts, frac)
  }
}

mask <- (fracts < 100) # filter out all values larger than 100 (could actually be already < N but so it shows that no other values lower then 100 were found)
hist(fracts[mask], breaks = seq(min(fracts[mask]) - 0.5, max(fracts[mask]) + 0.5, by = 1), xlab = "q", main = "Distribution of found q values smaller than 100") 

values <- sort(unique(fracts[mask])) # sort and remove duplicates

i = 1
while(values[i] <= N){
  r <- values[i]
  if ((3^r %% N == 1) && (r > 0)){ # check which q fulfills the equation 3^r mod N = 1 skip the 0 value, since it trivially fulfills the equation
    # stop the search if a solution is found
    print(r) # print the found solution 
    break
  }
  i <- i + 1
}

  



