library(qsimulatR)

x = qstate(3, coef = c(1, 2, 3, 4, 5, 6, 7, 8))
x = qstate(2)

print("initial state:")
print(x)

# x <- H(2)*x
x <- H(1)*x
x <- CNOT(c(1, 2))*x
print("after H(1):")
print(x)
