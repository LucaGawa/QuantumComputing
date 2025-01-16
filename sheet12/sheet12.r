library("qsimulatR")
library(hash)

create.state <- function() {
basis <- c()
for (i in c(0:(2^13-1))) {
  basis[i + 1] <-
    paste0("|helper=", 
           (i %/% 2^13) %% 2,
           (i %/% 2^12) %% 2,
           (i %/% 2^11) %% 2,
           (i %/% 1024) %% 2,
           (i %/% 256) %% 2,
           (i %/% 128) %% 2,
           ">|main=", 
           (i %/% 64) %% 2,
           (i %/% 32) %% 2,
           (i %/% 16) %% 2,
           (i %/% 8) %% 2, 
           (i %/% 4) %% 2,
           (i%/%2) %% 2, 
           i%%2, ">")
}
q <- qstate(13, basis = basis)
return(q)
}

encode <- function(q, bits) {
# creates the eigenstates of the stabilizers
# if the input ist |0000000> the output is |0_L>
# if the input ist |0000001> the output is |1_L>
# bits are 7 main bits of the state 
# in addition there should be 6 helper bits but they are theoretically not needed for the encoding

# map |0000001> into |1110000> and act as identity for |0000000>
q <- CNOT(c(bits[1], bits[7])) * q
q <- CNOT(c(bits[1], bits[6])) * q
q <- CNOT(c(bits[1], bits[5])) * q
q <- CNOT(c(bits[7], bits[1])) * q

# apply the unitary transformations

# (1 + g2g1)
q <- (CNOT(c(bits[3], bits[4])) * (CNOT(c(bits[3], bits[5])) * (CNOT(c(bits[3], bits[6])) * (H(bits[3]) * q))))
# (1 + g3g1)
q <- (CNOT(c(bits[2], bits[4])) * (CNOT(c(bits[2], bits[5])) * (CNOT(c(bits[2], bits[7])) * (H(bits[2]) * q))))
# q <- (CNOT(c(2,4)) * (CNOT(c(2,5)) * (CNOT(c(2,7)) * (H(2) * q))))
# (1 + g3g2g1)
q <- (CNOT(c(bits[1], bits[4])) * (CNOT(c(bits[1], bits[6])) * (CNOT(c(bits[1], bits[7])) * (H(bits[1]) * q))))
# q <- (CNOT(c(1,4)) * (CNOT(c(1,6)) * (CNOT(c(1,7)) * (H(1) * q))))

  return(q)
}

g1 <- function(x,bits) {
    x <- X(bits[1]) * x
    x <- X(bits[2]) * x
    x <- X(bits[3]) * x
    return(x)
  }

g2 <- function(x,bits) {
    x <- X(bits[1]) * x
    x <- X(bits[2]) * x
    x <- X(bits[5]) * x
    x <- X(bits[6]) * x
    return(x)
  }

g3 <- function(x,bits) {
    x <- X(bits[1]) * x
    x <- X(bits[3]) * x
    x <- X(bits[5]) * x
    x <- X(bits[7]) * x
    return(x)
  }

g4 <- function(x,bits) {
    x <- Z(bits[1]) * x
    x <- Z(bits[2]) * x
    x <- Z(bits[3]) * x
    x <- Z(bits[4]) * x
    return(x)
}

g5 <- function(x,bits) {
    x <- Z(bits[1]) * x
    x <- Z(bits[2]) * x
    x <- Z(bits[5]) * x
    x <- Z(bits[6]) * x
    return(x)
}

g6 <- function(x,bits) {
    x <- Z(bits[1]) * x
    x <- Z(bits[3]) * x
    x <- Z(bits[5]) * x
    x <- Z(bits[7]) * x
    return(x)
}


project.g1 <- function(q, main_bits, helper){
  # apply g1
  q <- g1(q, main_bits) 

  # change to correct basis
  q <- H(main_bits[1]) * q
  q <- H(main_bits[2]) * q
  q <- H(main_bits[3]) * q
  q <- H(main_bits[4]) * q
  
  # project to helper bit
  q <- CNOT(c(main_bits[1], helper)) * q
  q <- CNOT(c(main_bits[2], helper)) * q
  q <- CNOT(c(main_bits[3], helper)) * q
  q <- CNOT(c(main_bits[4], helper)) * q
  # q <- CNOT(c(main_bits[5], helper)) * q
  
  # change back to computational basis
  q <- H(main_bits[1]) * q
  q <- H(main_bits[2]) * q
  q <- H(main_bits[3]) * q
  q <- H(main_bits[4]) * q

  return(q)
}

project.g2 <- function(q, main_bits, helper){
  # apply g2
  q <- g2(q, main_bits) 

  # change to correct basis
  q <- H(main_bits[1]) * q
  q <- H(main_bits[2]) * q
  q <- H(main_bits[5]) * q
  q <- H(main_bits[6]) * q

  # project to helper bit
  q <- CNOT(c(main_bits[1],helper)) * q
  q <- CNOT(c(main_bits[2],helper)) * q
  q <- CNOT(c(main_bits[5],helper)) * q
  q <- CNOT(c(main_bits[6],helper)) * q
  
  # change back to computational basis
  q <- H(main_bits[1]) * q
  q <- H(main_bits[2]) * q
  q <- H(main_bits[5]) * q
  q <- H(main_bits[6]) * q

  return(q)
}

project.g3 <- function(q, main_bits, helper){
  # apply g3
  q <- g3(q, main_bits) 

  # change to correct basis
  q <- H(main_bits[1]) * q
  q <- H(main_bits[3]) * q
  q <- H(main_bits[5]) * q
  q <- H(main_bits[7]) * q

  # project to helper bit
  q <- CNOT(c(main_bits[1],helper)) * q
  q <- CNOT(c(main_bits[3],helper)) * q
  q <- CNOT(c(main_bits[5],helper)) * q
  q <- CNOT(c(main_bits[7],helper)) * q
  
  # change back to computational basis
  q <- H(main_bits[1]) * q
  q <- H(main_bits[3]) * q
  q <- H(main_bits[5]) * q
  q <- H(main_bits[7]) * q

  return(q)
}

project.g4 <- function(q, main_bits, helper){
  # apply g4
  q <- g4(q, main_bits) 

  # project to helper bit
  q <- CNOT(c(main_bits[1],helper)) * q
  q <- CNOT(c(main_bits[2],helper)) * q
  q <- CNOT(c(main_bits[3],helper)) * q
  q <- CNOT(c(main_bits[4],helper)) * q
  
  return(q)
}

project.g5 <- function(q, main_bits, helper){
  # apply g5
  q <- g5(q, main_bits) 

  # project to helper bit
  q <- CNOT(c(main_bits[1],helper)) * q
  q <- CNOT(c(main_bits[2],helper)) * q
  q <- CNOT(c(main_bits[5],helper)) * q
  q <- CNOT(c(main_bits[6],helper)) * q
  
  return(q)
}

project.g6 <- function(q, main_bits, helper){
  # apply g6
  q <- g6(q, main_bits) 

  # project to helper bit
  q <- CNOT(c(main_bits[1],helper)) * q
  q <- CNOT(c(main_bits[3],helper)) * q
  q <- CNOT(c(main_bits[5],helper)) * q
  q <- CNOT(c(main_bits[7],helper)) * q
  
  return(q)
}


gen.code.string <- function(q, main_bits) {

  # projects all states in the correct basis to the helper bits
  q <- project.g4(q, main_bits, 11)
  q <- project.g1(q, main_bits, 8)
  q <- project.g2(q, main_bits, 9)
  q <- project.g3(q, main_bits, 10)
  q <- project.g5(q, main_bits, 12)
  q <- project.g6(q, main_bits, 13)

  # measure the helper bits
  m1 <- measure(q, 8, rep=1)
  m2 <- measure(q, 9, rep=1)
  m3 <- measure(q, 10, rep=1)
  m4 <- measure(q, 11, rep=1)
  m5 <- measure(q, 12, rep=1)
  m6 <- measure(q, 13, rep=1)

  values <- c(m1$value, m2$value, m3$value, m4$value, m5$value, m6$value)

  string <- paste(sapply(values, function(mi) ifelse(mi == 1, "-1", "+1")), collapse = "") # translate form measurement result to eigenvalue (in theory everything could work just with 0 and 1 if used consistently, but since I wrote down the eigenvalues on the sheet in 2 and 3 I have choose this representation for a better comparison)

return(string)
}

find.eigenvalues <- function(plot=FALSE) {
  # calculates the all correction codes by measuring the helper bits systematically for all helper bits
  # if plot is TRUE the eigenvalue strings are printed, otherwise they get just stored in a hash table for later use

  main_bits <- c(1:7)
  h <- hash()
  h[["+1+1+1+1+1+1"]] <- Id(1)
  gates <- c(X, Y, Z)
  for (gate in gates) {
    for (i in main_bits) {
      q0 <- create.state()

      q <- encode(q0, main_bits)
      q <- gate(i) * q

      string <- gen.code.string(q, main_bits)
      h[[string]] <- gate(i)

      if (plot) {
        string <- paste0(gate(i)@type,"(", i, ") = ", string)
        print(string)
      }
    }
  }
  return(h)
}



h = find.eigenvalues(plot=TRUE) # precalculate all possible eigenvalues for the stabilizers and X,Y,Z Errors

########################
######## task 4 ########
########################

main_bits <- c(1:7)
q0 <- create.state() # create an 7+6 qubit state
#q0 <- X(1) * q0 # this one can be commented in in case one wants to create the state |1_L> instead of |0_L>

q <- encode(q0, main_bits) # encode the state in |0_L>
print("encoded")
print(q)

errors <- c("X", "Y", "Z") 
random_error <- sample(errors, 1) 
q <- noise(main_bits, error = random_error, type=random_error) * q # apply random an X,Y or Z error

print("error")
print(q)

string <- gen.code.string(q, main_bits) # find the eigenvalue string for the error
print(string)

q <- h[[string]] * q # correct the error according to the found eigenvalues 

print("corrected")
print(q)  

