library(qsimulatR)

print("Task 4:")
x4 = qstate(2)
print("Initial state:")
print(x4)
x4 <- H(2)*x4
x4 <- X(1)*x4
print("After applying H_2 and X_1:")
print(x4)

print("Task 5:")
x5 = qstate(2)
print("Initial state:")
print(x5)
x5 <- H(1)*x5
x5 <- CNOT(c(1, 2))*x5
print("After applying H_1 and CNOT_(1,2):")
print(x5)

print("More complex example:")
x = qstate(3, coef=c(1,2,3,4,5,6,7,8))
print("Initial state:")
print(x)
x <- H(1)*x
x <- CNOT(c(1, 2))*x
x <- CNOT(c(2, 3))*x
x <- H(3)*x
print("After applying H_1, CNOT_(1,2), CNOT_(2,3) and H_3:")
print(x)

