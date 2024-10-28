import numpy as np
import copy

class quant:
    def __init__(self, nqubits=2, coeffs=None, normalize=True):
        """ initialize a quantum state in the computational 
        basis. If no coefficients are given, the state 
        is set to |00...0> """

        self.n = nqubits
        self.state = np.zeros(2**self.n)
        
        if coeffs is not None:
            if len(coeffs) != 2**self.n:
                raise ValueError("Number of coefficients must be 2^n")
            if normalize:
                #normalize state to <psi|psi> = 1
                coeffs = coeffs / np.linalg.norm(coeffs)
            self.state[:len(coeffs)] = coeffs
        else:
            self.state[0] = 1
        pass
    
    def print_state(self, precision=3):
        """ print the quantum state """
        for i,coeff in enumerate(self.state):
            print(f"{coeff:.{precision}}|{bin(i)[2:]:>0{self.n}}>")
        return
    
    # gates
    def hadamard(self, qubit):
        """ apply Hadamard gate to qubit """
        pass

    def X(self, qubit):
        """ apply X (NOT) gate to qubit """
        old_state = copy.deepcopy(self.state)
        qubit_index = qubit - 1 # since q1 is the first qubit and not q0

        for i,coeff in enumerate(old_state):

            if (i >> qubit_index) & 1 == 0:
            # if (bin(i)[-qubit]) == '0':

                flipped_index = i | (1 << qubit_index)

                self.state[i] = old_state[flipped_index]
                self.state[flipped_index] = old_state[i]
        print(old_state)


if __name__ == "__main__":
    coeffs = np.arange(1,9)
    print(coeffs)
    x = quant(3, coeffs=coeffs)
    
    print("Initial state:")
    x.print_state(precision=5)
    print("apply X gate to qubit 1:")
    x.X(2)
    x.print_state(precision=5)
    
