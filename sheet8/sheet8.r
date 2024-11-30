library('qsimulatR')

n <- 6
psi <- qstate(n)
psi <- X(n) * psi # apply a delta peak in position space
# for (i in 1:(n)) {
#   psi <- X(i) * psi
#   psi <- H(i) * psi
# }
# psi <- H(2) * psi
# mes = measure(psi) # measure the state
# plot(hist(mes)) # plot the histogram of the state

Rtheta <- function(i, theta) {
  return (sqgate(bit=i, M=array(as.complex(c(1,0,0,exp(1i*theta))), dim=c(2,2))))
}


free.time.evolution <- function(psi, dt) {
n <- psi@nbits
dp <- 2*pi/(2^n)
a <- pi^2/(n^2)
b <- pi*dp/n
c <- dp^2

q <- qft(psi)

# apply the global phase
theta <- -n^2*a*dt/2
q <-  Rtheta(1,theta)*(X(1)*(Rtheta(1, theta) * (X(1) * q)))

for (j in 1:n) {
  theta <- dt/2*(n*b*2^(j)-c*2^(2*j-2))

  q <- Rtheta(j, theta) * q
  for (l in 1:n){
    if (l != j) {
      theta <- -c*dt/2*2^(l+j-2)
      q <- cqgate(bits=c(l,j), gate=Rtheta(j, theta)) * q
    }
  }
}
q <- qft(q, inverse=TRUE)
return(q)
}

harmonic.potential <- function(psi, dt, gamma, L){
  n <- psi@nbits
  dx <- L/(2^n)
  a <- gamma*(L/2)^2/(n^2)
  b <- gamma*L/2*dx/n
  c <- gamma*dx^2

# apply the global phase
theta <- -n^2*a*dt/2
q <-  Rtheta(1,theta)*(X(1)*(Rtheta(1, theta) * (X(1) * psi)))

for (j in 1:n) {
  theta <- dt/2*(n*b*2^(j)-c*2^(2*j-2))

  q <- Rtheta(j, theta) * q
  for (l in 1:n){
    if (l != j) {
      theta <- -c*dt/2*2^(l+j-2)
      q <- cqgate(bits=c(l,j), gate=Rtheta(j, theta)) * q
    }
  }
}
  return(q)
}

# mes <- measure(psi, rep=1000)
# plot(hist(mes,only.nonzero=FALSE, freq=FALSE))

dt <- 1
for (i in (1:10)) {
  psi <- harmonic.potential(psi, dt, 1, 1000)
  psi <- free.time.evolution(psi, dt)
  if (i %% 1 == 0) {
  mes <- measure(psi, rep=100000)
  plot(hist(mes,only.nonzero=FALSE, freq=FALSE))
  }
}



# print(psi)
