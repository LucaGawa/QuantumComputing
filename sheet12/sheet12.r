library("qsimulatR")

# basis <- c()
n <- 2^13
basis <- character(n)
for(i in c(0:(2^13-1))) {
  basis[i + 1] <-
    paste0("|helper=", sprintf("%06d", as.integer(intToBits(i %/% 2^7)[1:6])),
           ">|main=", sprintf("%07d", as.integer(intToBits(i %% 2^7)[1:7])),
           ">")
}


encode <- function(q, bits) {
q <- CNOT(c(bits[1], bits[7])) * q
q <- CNOT(c(bits[1], bits[6])) * q
q <- CNOT(c(bits[1], bits[5])) * q
q <- CNOT(c(bits[7], bits[1])) * q


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
# map |0000001> into |1110000>

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

Sdagger <- function(bit) {
  # implements the dagger of the S gate
  methods::new("sqgate", bit=as.integer(bit), M=array(as.complex(c(1, 0, 0, -1i)), dim=c(2,2)),
  type="Sd")
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

#
# # q0 <- qstate(13, basis = basis)
# # q1 <- qstate(13, basis = basis)
q0 <- qstate(13)
q1 <- qstate(13)
q1 <- X(1) * q1
main_bits <- c(1:7)

print("q0")
print(q0)
print("q1")
print(q1)

q <- encode(q0, c(1:7))

print("encoded")
print(q)

plot.eigs <- function(gate) {
for (i in main_bits) {
  q0 <- qstate(13)
  q <- encode(q0, main_bits)
  q <- gate(i) * q

  q <- project.g4(q, main_bits, 11) #todo check why order is relevant
  q <- project.g1(q, main_bits, 8)
  q <- project.g2(q, main_bits, 9)
  q <- project.g3(q, main_bits, 10)
# q <- project.g4(q, main_bits, 11)
  q <- project.g5(q, main_bits, 12)
  q <- project.g6(q, main_bits, 13)


m1 <- measure(q, 8, rep=1)
m2 <- measure(q, 9, rep=1)
m3 <- measure(q, 10, rep=1)
m4 <- measure(q, 11, rep=1)
m5 <- measure(q, 12, rep=1)
m6 <- measure(q, 13, rep=1)

values <- c(m1$value, m2$value, m3$value, m4$value, m5$value, m6$value)
# values <- c(m1, m2, m3, m4, m5, m6)

string <- paste(sapply(values, function(mi) ifelse(mi == 1, "-1", "+1")), collapse = "")

string <- paste0(as.character(substitute(gate)),"(", i, ") = ", string)
print(string)

}
}


plot.eigs(X)
plot.eigs(Z)
plot.eigs(Y)

# print("projected")
# print(q)
# q <- g1(q, c(1:7))
# plot(q)
# print("")

# print(q)
# # print()  
# # print(
# # summary(res)
# # q <- noise(1:7, error = "X", type="X") * q
# #
# # # correct.flip <- function(q, kj
# #
# # q <- qstate(21)
# # print(q)
# #
# # plot(q)
