library('qsimulatR')

#addition by fourier transform

n <- 3 

x = 2

coefs <- rep(0, 2^n)
# coefs[58] <- 1
coefs[4] <- 1
# coefs <- rep(0, 2^(2*n))
# coefs[2] = 1

q <- qstate(n, coefs=coefs, basis=as.character(0:(2^n-1)))

print(q)
print('')

q <- qft(q)


for (j in 1:n) {

    theta = 2 * pi *  x / (2^(n-j+1))
    q <- sqgate(bit=j, M=array(as.complex(c(1,0,0,exp(1i*theta) )), dim=c(2,2))) * q

}

q <- qft(q, inverse=TRUE)
print(q)
