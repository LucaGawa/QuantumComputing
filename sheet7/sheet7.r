library('qsimulatR')

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

addmodN <- function(psi, y, N, a, c1, c2, xbits){
  
  n <- length(xbits)
  
  # make sure that y is smaller than N
  y <- y %% N

  q <- is.less( psi, y=N, a=a, c1=c1, xbits=xbits) #check if x < N and store the result in c1
  q <- is.less( q, y=N-y, a=a, c1=c2, xbits=xbits) #check if x < N-y and store the result in c2


    
  # if c1=1 and c2=0 compute x=x+y-N
  q <- X(c2) * (CCNOT(c(c1, c2, a))* (X(c2)*q))
  q <- cadd(a, q, y=y-N, xbits=xbits) 
  q <- X(c2)* (CCNOT(c(c1, c2, a))*(X(c2)*q))


  # if c1=1 and c2=1 compute x=x+y
  q <- CCNOT(c(c1, c2, a))*q
  q <- cadd(a, q, y=y, xbits=xbits)
  q <- CCNOT(c(c1, c2, a))*q
  
  #reset helper bits
  q <- is.less(q, y=y, a=a, c1=c2, xbits=xbits)
  q <- CNOT(c(c1, c2)) * q
  q <- is.less(q, y=N, a=a, c1=c1, xbits=xbits)
  
  return(q)
}

creat.qstate <- function(nbits){
  # create a quantum state with nbits which represent the register x with length nbits-3 and
  # 3 additional helper bits a, c2, c1, which are needed for the addmodN function
  # |x>|a>|c2>|c1>
  # returns a qstate object initialized in the state |0>|0>|0>|0> 
  basis <- c()
  for (i in 0:(2^nbits - 1)) {

  binary <- intToBits(i)[1:nbits] # Convert to binary 
  binary <- rev(as.integer(binary)) # fit the order to the convention
  
  left_bits <- sum(binary[1:(nbits - 3)] * 2^(rev(0:(nbits - 4))))
  right_bits <- binary[(nbits - 3):nbits]
  
  basis[i + 1] <- paste0(
    "|", left_bits, ">|a=", right_bits[1], 
    ">|c2=", right_bits[2], ">|c1=", right_bits[3], ">"
  )
}

psi <- qstate(nbits, basis=basis) 
return(psi)
}


print("#### Task 3 ####")
x = 7
y = 2
N = 10
nbits <- 10 # has to be large enough such that N <= 2^n. In additional there has to enough bits the 3 helper bits
psi <- creat.qstate(nbits) # initialize the quantum state with the necessary bits 
psi <- add(psi, x, xbits=c(4:nbits)) # add a value to the register x to make it unequal to 0
print("state before addmodN")
print(psi)
psi <- addmodN(psi, y=y, N=N, a=3, c1=1, c2=2, xbits=c(4:nbits))
print("state after addmodN")
print(psi)
print(sprintf("theoretical result (without overflow): %d", (x+y) %% N))

