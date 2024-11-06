library(qsimulatR)

################
#    Task 6    #
################

# initial state |0>x|0^3>
x <- qstate(4)

# apply H on |0>
x <- X(4)*x

# apply H on every qubit
for (i in 1:4) {
  x <- H(i)*x
}

# apply Us for s=010 -> CNOT(2,4)
x <- CNOT(c(2,4))*x

# apply H on first n=3 qubits
for (i in 1:3) {
  x <- H(i)*x
}

plot(x)
# print(x)

################
#    Task 7    #
################

# create state |psi(x=1/4)>

# starting with state |0000>
q <- qstate(4)

# apply H on every qubit
for (i in 1:4) {
  q <- H(i)*q
}

# apply Rtheta with the correct theta to every qubit
for (i in 1:4) {
  q <- Ri(i, 7-i, sign=-1)*q
}

################
#    Task 8    #
################

psi <- qft(q)
rep <- 1000000000

hist(measure(psi, rep=rep), freq=FALSE)
summary(measure(psi, rep=rep))


