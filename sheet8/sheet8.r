library('qsimulatR')


# define the gate R(theta) for applying a phase to a single qubit
Rtheta <- function(i, theta) {
  return (sqgate(bit=i, M=array(as.complex(c(1,0,0,exp(1i*theta))), dim=c(2,2))))
}


free.time.evolution <- function(psi, dt) {
  # perform the free time evolution for a single trotter time step dt
  # input psi is the quantum state
  

n <- psi@nbits
dp <- 2*pi/(2^n) # discretization of the momentum space
a <- pi^2/(n^2)
b <- pi*dp/n
c <- dp^2

# transform to the momentum space
q <- qft(psi) # for a better performance and without a potential one do not need to transform the state for every time step back and forth (for the compacter code in the loop this is not done here for this excercise)

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

# transform back to the position space
q <- qft(q, inverse=TRUE)
return(q)
}

harmonic.potential <- function(psi, dt, gamma, L){
  # perform the harmonic potential part of the time evolution for a single trotter time step dt
  # input psi is the quantum state

  n <- psi@nbits
  dx <- L/(2^n) # discretization of the position space
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


n <- 6
psi <- qstate(n)
psi <- X(n) * psi # apply a delta peak in position space

dt <- 1 # for a good approximation this should be smaller, but for showing the effect in a rough simulation this seems to be enough

# task 2.1
for (i in (1:11)) {
  mes <- measure(psi, rep=100000)
  hist(mes,only.nonzero=FALSE, freq=FALSE, ylim=c(0,1))
  psi <- free.time.evolution(psi, dt) # apply the kinetic part in the momentum space
}

# task 2.3
psi <- qstate(n)
psi <- X(n) * psi # apply a delta peak in position space
for (i in (1:10)) {
  mes <- measure(psi, rep=100000)
  hist(mes,only.nonzero=FALSE, freq=FALSE, ylim=c(0,1))
  psi <- free.time.evolution(psi, dt) # apply the kinetic part in the momentum space
  psi <- harmonic.potential(psi, dt, 100, 10) # apply the potential part in the real space
}

