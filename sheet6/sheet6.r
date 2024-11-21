library('qsimulatR')

add_by_fourier <- function(x,y) {
  #Addition of the two registers x and y. 
  #The input of x and y are strings in the form e.g. '011' or '1110' both must have
  #the same length and get converted to a quantum state |x,y>.
  n <- nchar(x) 
  ntest <- nchar(y)
  if (n != ntest) {
    stop("register x and register y have to have the same length n")
  }
  #initialize a quantum state with |x,y>
  state = paste(x, y, sep='')
  decimal_value <- strtoi(state, base = 2)

  coefs <- rep(0, 2^(2*n))
  coefs[decimal_value+1] = 1

  q <- qstate(2*n, coefs=coefs)
  
  print("the input state |x,y> is:")
  print(q)
  decimal_y <- strtoi(y, base = 2)
  decimal_x <- strtoi(x, base = 2)
  formatted_string <- sprintf("This corresponds to (%d + %d) mod 2^%d = %d", decimal_x, 
                              decimal_y,n, (decimal_y + decimal_x) %% 2^n)
  print(formatted_string)

  # main part
  
  # apply fourier transform on y
  q <- qft(q, bits=c(1:n))

  # apply controlled phase gates
  for (j in 1:n) {
    for (k in (n+1):(2*n)) {
    q <- cqgate(bits=c(k, j), gate=Ri(bit=j, i=(-j-k+2*n+2))) * q
  }
}
  
  # apply inverse fourier transform on y
  q <- qft(q, bits=c(1:n), inverse=TRUE)
  print("the output state |x,y> is:")
  print(q)
  # plot(q)
}


# addition and mulitplication by fourier transform
add_mul_by_fourier <- function(x,y,a) {
  #addition and multiplication by fourier transform: |x,y> -> |x+ay, y mod 2^n>
  #The input of x and y are strings in the form e.g. '011' or '1110' both must have
  #the same length and get converted to a quantum state |x,y>. a is an integer.
  n <- nchar(x) 
  ntest <- nchar(y)
  if (n != ntest) {
    stop("register x and register y have to have the same length n")
  }
  #initialize a quantum state with |x,y>
  state = paste(x, y, sep='')
  decimal_value <- strtoi(state, base = 2)

  coefs <- rep(0, 2^(2*n))
  coefs[decimal_value+1] = 1

  q <- qstate(2*n, coefs=coefs)
  
  print("the input state |x,y> is:")
  print(q)
  decimal_y <- strtoi(y, base = 2)
  decimal_x <- strtoi(x, base = 2)
  formatted_string <- sprintf("This corresponds to (%d + %d*%d) mod 2^%d = %d", decimal_x, a, 
                              decimal_y,n, (decimal_x + a*decimal_y) %% 2^n)
  print(formatted_string)

  # main part
  
  # apply fourier transform on x
  q <- qft(q, bits=c((n+1):(2*n)))

  # apply controlled phase gates
  for (j in (n+1):(2*n)) {
    for (k in 1:n) {

    theta = 2 * pi * a  / (2^(2*n-j-k+2)) # this is the same phase then before just withe the integer a added
    phasegate <- sqgate(bit=j, M=array(as.complex(c(1,0,0,exp(1i*theta) )), dim=c(2,2))) 
    q <- cqgate(bits=c(k, j), gate=phasegate) * q
  }
}
  
  # apply inverse fourier transform on x
  q <- qft(q, bits=c((n+1):(2*n)), inverse=TRUE)
  print("the output state |x,y> is:")
  print(q)
  # plot(q)
}


print("#### Task 4 ####")
add_by_fourier('111', '011')
print("################")
print("#### Task 6 ####")
add_mul_by_fourier('0010', '0011', 3)
print("################")

