library("qsimulatR")

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




# cmult <- function(c, psi, y, xbits){
#   n <- length(xbits)
#   q <- cqft(c, psi, bits=xbits)

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

# caddmodN <- function(c, psi, y, N, a, c1, c2, xbits){

cmultmodN <- function(c, psi, y, N, xbits, reg2, helpers){
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

# cexpmodN <- function(c, j, x, y, N, xbits, reg2, helpers){
#   # controlled version of |x>|0> -> |x^j*y mod N>|0>
#   # xbits has to have the same length as reg1
#   # 4 helper bits are required
#   psi <- x
#   for (i in c(1:j)) {
#     psi <- cmultmodN(c, psi, y, N, xbits, reg2, helpers)
#   }
#   return(psi)
# }

cexpmodN <- function(c, j, x, y, N, xbits, reg2, helpers){
  # controlled version of |x>|0> -> |x^j mod N>|0>
  # xbits has to have the same length as reg1
  # 4 helper bits are required
  psi <- x
  ab <- as.integer(intToBits(j))
  n <- max(which(ab == 1))
  y2 <- y %% N
  for (i in c(1:n)) {
    if (ab[i] == 1) {
    psi <- cmultmodN(c, psi, y=y2, N=N, xbits=xbits, reg2=reg2, helpers=helpers)
    }
    y2 <- ((y2 %% N) * (y2 %% N)) %% N
  }
  return(psi)
}



generate_basis <- function(n) {
  basis <- c()
  total_states <- 2^(2 * n + 5 + (2 * n + 3))  # Total states include `t` register (2n+3 qubits), `reg1`, `reg2`, helpers, and control bit
  
  for (i in 0:(total_states - 1)) {
    t <- i %/% (2^(2 * n + 5))  # Extract `t` (most significant qubits)
    xbits <- (i %/% (2^(n + 5))) %% 2^n  # Extract `reg1`
    reg2 <- (i %/% (2^5)) %% 2^n  # Extract `reg2`
    helpers <- c((i %/% 16) %% 2, (i %/% 8) %% 2, (i %/% 4) %% 2, (i %/% 2) %% 2)  # Extract the 4 helper qubits
    c_bit <- i %% 2  # Extract the control bit
    
    basis[i + 1] <- paste0(
      "|t=", t, 
      ">|xbits=", xbits, 
      ">|reg2=", reg2, 
      ">|helpers=", paste(helpers, collapse = ""), 
      ">|c=", c_bit, ">"
    )
  }
  
  return(basis)
}
n <- 3
basis <- generate_basis(n)
bits_total <- 4*n+7+1 
tbits_total <- 2*n+3

# q <- qstate(bits_total, basis=basis)
# q <- X(1) * q
# q <- X(9) * q
# # q <- X(9) * q
#
# print(q)

# q <- phase_estimation(bitmas=c(bits_total-tbits_total:bits_total),
#                       FUN=cexpmodN, x=q, y=1, N=7, xbits=bits_total-tbits_total-n:bits_total, reg2=bits_total-tbits_total-2*n:bits_total-tbits_total-n, helpers=bits_total-tbits_total-2*n-4:bits_total-tbits_total-2*n)
# # q <- X(6) * q
# print(q)

# q <- H(1) * q
# q <- CNOT(c(1,10)) * q

# ##### 
# # n <- 2
# #####
# c <- 1
# helpers <- c(2:5)
# reg2 <- c(6:7)
# xbits <- c(8:9)
# bitmas <- c(10:bits_total)
# N <- 3
# y <- 2

#####
n <- 3
#####
c <- 1
helpers <- c(2:5)
reg2 <- c(6:8)
xbits <- c(9:11)
N <- 7
y <- 3
bitmas <- c(12:bits_total)

# q <- phase_estimation(bitmas=bitmas,
#                       FUN=cexpmodN, x=q, y=y, N=N, xbits=xbits, reg2=reg2, helpers=helpers)

# phi = 0
# for (bit in bitmas) {
#   res <- summarmeasure(q, repetitions=5000, bit=bit)
#   print(res)
#   # resbits <- genStateNumber(whic
#   # print(resbits)
# }


# res = summary(measure(q, repetitions=5000))
# resbits <- genStateNumber(which(res$value==1)-1, tbits_total+1)
# phi <- sum(resbits[1:tbits_total] / 2^(1:tbits_total))
# print(phi)

# print(q)
# hist(res)
# res <- cmultmodN(c, q, y=2, N=N, xbits=xbits, reg2=reg2, helpers=helpers)
# print(res)

# res <- cmultmodN(c, res, y=3, N=N, xbits=xbits, reg2=reg2, helpers=helpers)
# print(res)
# res <- cmultmodN(c, res, y=3, N=N, xbits=xbits, reg2=reg2, helpers=helpers)
# print(res)
# res <- cmultmodN(c, res, y=3, N=N, xbits=xbits, reg2=reg2, helpers=helpers)
# print(res)

# q <- X(1) * q
# q <- caddmodN(c=1, psi=q, y=5, N=12, a=4, c1=2, c2=3, xbits=5:n)
# q <- caddmodN(c=1, psi=q, y=4, N=7, a=4, c1=2, c2=3, xbits=5:n)
# print(q)




# print(summands(3, 2^3, 3))

continued_fraction <- function(x, eps=1e-14, k_max=100){
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

     # Break condition

  
      if (abs(x0 - sn) < eps) {
        break
      }

    x <- 1/frac_part
}
  print(q)
  return(k)
}

print(continued_fraction(79/2^9))



