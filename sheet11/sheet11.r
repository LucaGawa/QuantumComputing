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

create_basis <- function(n) {
  basis <- c()
  for(i in 0:(2^n - 1)) {
    # binary <- intToBits(i)[1:n]
    basis[i + 1] <- 
      paste0("|", i %/% 16, ">|a=", (i %/% 8) %% 2, 
             ">|c2=", (i %/% 4) %% 2, ">|c1=", (i %/% 2) %% 2, ">|c=", i %% 2, ">")
  }
  return(basis)
}

n <- 10
basis <- create_basis(n)
q <- qstate(n, basis=basis)

q <- X(1) * q
q <- caddmodN(c=1, psi=q, y=5, N=12, a=4, c1=2, c2=3, xbits=5:n)
q <- caddmodN(c=1, psi=q, y=4, N=7, a=4, c1=2, c2=3, xbits=5:n)
print(q)






