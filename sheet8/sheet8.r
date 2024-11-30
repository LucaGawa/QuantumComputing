library('qsimulatR')

n <- 6
psi <- qstate(n)
psi <- X(n) * psi # apply a delta peak in position space
# mes = measure(psi) # measure the state
# plot(hist(mes)) # plot the histogram of the state

Rtheta <- function(i, theta) {
  return (sqgate(bit=i, M=array(as.complex(c(1,0,0,exp(1i*theta))), dim=c(2,2))))
}

# psi <- SWAP(c(1,n)) * psi
# psi <- SWAP(c(2,n-1)) * psi
# psi <- SWAP(c(3,n-2)) * psi

time.evolution <- function(psi, dt) {
# dt <- 0.1
n <- psi@nbits
dp <- 2*pi/(2^n)
a <- pi^2/(n^2)
b <- pi*dp/n
c <- dp^2
# a <- pi^2
# b <- pi^2/2^(n-1)
# c <- (pi^2/2^n)^n

q <- qft(psi)

theta <- -n^2*a*dt/2
# theta <- -dt/2*a
#
q <-  Rtheta(1,theta)*(X(1)*(Rtheta(1, theta) * (X(1) * q)))
# q <- X(1) * q
# q <- Rtheta(1, theta) * q
# q <- X(1) * q
# q <- Rtheta(1, theta) * q

for (j in 1:n) {
  # theta <- (2*b*2**(n-j-1)-c*2**(2*n-2-2*j))*dt/2
  # theta <- n*dt*b*2^(j-1)-dt/2*c*4^(j-1)
  theta <- 2*n*b*dt/2*2^(j-1)-c*dt/2*2^(2*j-2)

  q <- Rtheta(j, theta) * q
  for (l in 1:n){
    if (l != j) {
      # theta <- -dt*c*2^(l+j-2)
      # theta <- -c*2**(2*n-2-l-j)*dt/2
      theta <- -c*dt/2*2^(l+j-2)
      q <- cqgate(bits=c(l,j), gate=Rtheta(j, theta)) * q
    }
  }
}


q <- qft(q, inverse=TRUE)
  return(q)
}

# mes <- measure(psi, rep=1000)
# plot(hist(mes,only.nonzero=FALSE, freq=FALSE))

for (i in (1:100)) {
  psi <- time.evolution(psi, 0.1)
mes <- measure(psi, rep=1000)
plot(hist(mes,only.nonzero=FALSE, freq=FALSE))
}



# print(psi)
