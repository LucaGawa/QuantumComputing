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

    def X(self, qubit):
        """ apply X (NOT) gate to qubit """
        old_state = copy.deepcopy(self.state)

        for i,coeff in enumerate(old_state):

            # state = list(format(i, f'0{self.n}b')) # convert to binary 
            # flip_index = (int(state[-qubit])+1)%2 # 0 -> 1, 1 -> 0
            # state[-qubit] = str(flip_index) # apply the changed bit to the state
            # index = int("".join(state),2) # convert back to decimal

            index = i ^ (1 << (qubit - 1))  # finding the "partner index"

            # exchange the coefficients of the corresponding states
            self.state[i] = old_state[index] 
            self.state[index] = old_state[i]

        return

    def H(self, qubit):
        """ apply Hadamard gate to qubit """
        old_state = copy.deepcopy(self.state)

        for i,coeff in enumerate(old_state):
            
            index = i ^ (1 << (qubit - 1)) # 

            # apply the local operation
            self.state[i] = (old_state[index]-old_state[i]) / np.sqrt(2)
            self.state[index] = (old_state[i] + old_state[index]) / np.sqrt(2)
        return




if __name__ == "__main__":
    coeffs = np.arange(1,9)
    print(coeffs)
    x = quant(3, coeffs=coeffs)
    
    print("Initial state:")
    x.print_state(precision=5)
    print("apply X gate to qubit 1:")
    x.H(1)
    x.print_state(precision=5)
    
