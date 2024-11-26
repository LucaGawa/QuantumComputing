library('qsimulatR')

addmodN <- function(x, y, N) {
  # Add the integer y to the register x modulo the integer N
  # x has to be chosen large enough to represent the sum x+y mod N

  n <- nchar(x)
  # initialize a quantum state with |x>
  decimal_value <- strtoi(x, base = 2)

  coefs <- rep(0, 2^n)
  coefs[decimal_value + 1] <- 1
  q <- qstate(n, coefs=coefs)

  
  print("the input state |x> is:")
  print(q)

  q <- qft(q)

  for (j in 1:n) {
  theta <- 2 * pi * y * 2^(j-1-n)/ N 
  phasegate <- sqgate(bit=j, M=array(as.complex(c(1,0,0,exp(1i*theta) )), dim=c(2,2))) 
  # q <- cqgate(bits=c(k, j), gate=phasegate) * q
  q <- phasegate * q
  }

  q <- qft(q, inverse=TRUE)
  
  print("the output state |x+y mod N> is:")
  print(q)
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

cis.less <- function(c, phi, y, a, c1, xbits){
  
  bits <- c(xbits, a) # add the helper bit for the underflow
  n <- length(bits)

  q <- cadd(c, psi, y=2^n - y, xbits=bits) # subtract y 
  q <- CNOT(c(a, c1)) * q # store the value of the underflow to the helper qubit c1 

  q <- cadd(c, q, y, bits) # add y back to restore the original state beside the value of c1
  return(q)
}


basis <- c()
for(i in c(0:(2^6-1))) {
  basis[i + 1] <-
    paste0("|", i %/% 8 ,">|a=",
    (i %/% 4) %% 2, ">|c1=", (i%/%2) %% 2,
    ">|c=", i%%2, ">")
}

psi <- qstate(6,basis=basis) # |x>|a>|c1>|c>
psi <- X(1) * psi
psi <- cadd(1, psi, 3, xbits=c(4:5))
print(psi)
psi <- cis.less(1, x, 5, a=3, c1=2, xbits=c(4:5))
# psi <- cadd(1, psi, 3, xbits=c(2:4) )
# psi <- X(1) * psi

print(psi)

# print("#### Task 3 ####")
# addmodN('0101', 1, 2^4)
