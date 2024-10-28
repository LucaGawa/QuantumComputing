import numpy as np

class quant:
    def __init__(self,nqubits=2, coeffs=None, normalize=True):
        """ initialize a quantum state in the computational 
        basis. If no coefficients are given, the state 
        is set to |00...0> """

        self.n = nqubits
        self.state = np.zeros(2**self.n)
        
        if coeffs is not None:
            if len(coeffs) != 2**self.n:
                raise ValueError("Number of coefficients must be 2^n")
            if normalize:
                coeffs = coeffs / np.linalg.norm(coeffs)
            self.state[:len(coeffs)] = coeffs
        else:
            self.state[0] = 1
        pass
    
    def print_state(self):
        """ print the quantum state """
        for i,coeff in enumerate(self.state):
            print(f"{coeff:.3f}|{bin(i)[2:]:>0{self.n}}>")
        return


coeffs = np.full(2**3,1)
print(coeffs)
x = quant(3, coeffs=coeffs)

x.print_state()
    
