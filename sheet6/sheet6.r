library('qsimulatR')

#addition by fourier transform numer

n <- 8 

x = 2

coefs <- rep(0, 2^n)
# coefs[58] <- 1
coefs[4] <- 1
# coefs <- rep(0, 2^(2*n))
# coefs[2] = 1

q <- qstate(n, coefs=coefs, basis=as.character(0:(2^n-1)))

# print(q)
# print('')

q <- qft(q)


for (j in 1:n) {

    theta = 2 * pi *  x / (2^(n-j+1))
    q <- sqgate(bit=j, M=array(as.complex(c(1,0,0,exp(1i*theta) )), dim=c(2,2))) * q
}

q <- qft(q, inverse=TRUE)
# print(q)






########################################################3


# addition by fourier transform two registers
n <- 3

coefs <- rep(0, 2^(2*n))
# coefs[33] = 1
coefs[3] = 1
q <- qstate(2*n, coefs=coefs)
print(q)

q <- qft(q, bits=c(1:n))

for (j in 1:n) {
    # print(j)
    print("")
  for (k in 4:6) {
    # print(k)
    theta = 2 * pi *  2^(k-n-1)  / (2^(n-j+1))
    phasegate <- sqgate(bit=j, M=array(as.complex(c(1,0,0,exp(1i*theta) )), dim=c(2,2))) 
    q <- cqgate(bits=c(k, j), gate=phasegate) * q
    # print(exp(1i*theta))
  }
}

q <- qft(q, bits=c(1:n), inverse=TRUE)

plot(q)

print(q)


########################################################

# multiply and add 
a <- 3
x <- 2

coefs <- rep(0, 2^n)
coefs[4] <- 1

q <- qstate(n, coefs=coefs, basis=as.character(0:(2^n-1)))

# print(q)
# print('')
q <- qft(q)

for (j in 1:n) {

    theta = 2 * pi *  x / (2^(n-j+1))
    q <- sqgate(bit=j, M=array(as.complex(c(1,0,0,exp(1i*theta*a) )), dim=c(2,2))) * q
}

q <- qft(q, inverse=TRUE)
# print(q)



