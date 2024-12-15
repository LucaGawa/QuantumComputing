library("qsimulatR")

create.qstate <- function(n, representation = "decimal") {
  # Create a quantum state with n qubits representing the register x
  # and 1 additional helper qubit at the far right
  # |x>|h>
  # Returns a qstate object initialized in the state |0>|0>
  # 
  # Parameters:
  #   n: number of qubits in the main register
  #   representation: "decimal" or "binary" for the main register representation
  
  basis <- c()
  for (i in 0:(2^(n+1) - 1)) {
    binary <- intToBits(i)[1:(n+1)]  # Convert to binary
    binary <- rev(as.integer(binary))  # Fit the order to the convention
    
    if (representation == "decimal") {
      main_register <- sum(binary[1:n] * 2^(rev(0:(n-1))))
    } else if (representation == "binary") {
      main_register <- paste(binary[1:n], collapse="")
    } else {
      stop("Invalid representation. Use 'decimal' or 'binary'.")
    }
    
    helper_qubit <- binary[n+1]
    
    basis[i + 1] <- paste0("|", main_register, ">|", helper_qubit, ">")
  }
  
  psi <- qstate(n+1, basis=basis)
  return(psi)
}

Sdagger <- function(bit) {
  # implements the dagger of the S gate
  methods::new("sqgate", bit=as.integer(bit), M=array(as.complex(c(1, 0, 0, -1i)), dim=c(2,2)),
  type="Sd")
}

dt <- 1
apply.H <- function(psi, dt){
n <- 3
psi <- create.qstate(n, representation="binary")
print(psi)

# apply P transformations 
psi = S(2)*(H(2)*psi)
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
psi = H(2)*(Sdagger(2)*psi)
return(psi)
}


# apply measurement
psi = apply.H(psi, dt)
mes <- measure(psi, rep=10000)
hist(mes, only.nonzero=TRUE, freq=FALSE)

plot(psi)
print(psi)

