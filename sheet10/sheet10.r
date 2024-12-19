library("qsimulatR")

Sdagger <- function(bit) {
  # implements the dagger of the S gate
  methods::new("sqgate", bit=as.integer(bit), M=array(as.complex(c(1, 0, 0, -1i)), dim=c(2,2)),
  type="Sd")
}

dt <- 1
n <- 3
psi <- qstate(4)
apply.H <- function(psi, dt){

# apply P transformations 
psi = H(2)*(S(2)*psi)
psi = H(3)*psi
# apply CNOT gates
for (i in 1:3) {
  psi = CNOT(c(i,4))*psi
}
# apply phase gate
psi = Rz(4, theta=2*dt)*psi

# apply CNOT gates
for (i in 1:3) {
  psi = CNOT(c(n+1-i,4))*psi
}
# apply Pdagger transformations
psi = H(3)*psi
psi = Sdagger(2)*(H(2)*psi)
return(psi)
}

# apply measurement
psi = apply.H(psi, dt)
print(psi)
mes <- measure(psi, rep=10000)
hist(mes, only.nonzero=TRUE, freq=FALSE)
plot(psi)
