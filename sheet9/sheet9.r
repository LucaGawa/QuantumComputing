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
n = 3
q = create.qstate(n, representation="binary")
for (i in 2:(n+1)) {
  q <- H(i)*q
}
for (i in 2:(n+1)) {
  q <- X(i)*q
}
q <- cnqgate(cbits=c(2:(n+1)), tbit=1, gate=X(1)) * q
for (i in 2:(n+1)) {
  q <- X(i)*q
}
# print("after")
# print(q)
theta <- 1
phasegate <- sqgate(bit=1, M=array(as.complex(c(1,0,0,exp(1i*theta) )), dim=c(2,2))) 
q <- phasegate * q
for (i in 2:(n+1)) {
  q <- X(i)*q
}
print("")
print(q)
print("")
q <- cnqgate(cbits=c(2:(n+1)), tbit=1, gate=X(1)) * q
for (i in 2:(n+1)) {
  q <- X(i)*q
}

for (i in 2:(n+1)) {
  q <- H(i)*q
}
print(q)
