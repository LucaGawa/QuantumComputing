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
    
    def print_state(self, precision=3, print_zeros=False):
        """ print the quantum state """
        for i,coeff in enumerate(self.state):
            if not print_zeros and np.isclose(coeff, 0):
                continue
            print(f"{coeff:.{precision}}|{bin(i)[2:]:>0{self.n}}>")
        return
    
    # gates
    def X(self, qubit):
        """ apply X (NOT) gate to qubit """
        old_state = copy.deepcopy(self.state) # avoid overwriting

        for i,_ in enumerate(old_state):

            index = i ^ (1 << (qubit - 1))  # finding the "partner index"

            # exchange the coefficients of the corresponding states
            self.state[i] = old_state[index] 
            self.state[index] = old_state[i]

        return self

    def H(self, qubit):
        """ apply Hadamard gate to qubit """
        old_state = copy.deepcopy(self.state) # avoid overwriting

        for i,_ in enumerate(old_state):
            
            index = i ^ (1 << (qubit - 1)) # finding the "partner index"

            # apply the local operation
            self.state[i] = (old_state[index]-old_state[i]) / np.sqrt(2)
            self.state[index] = (old_state[i] + old_state[index]) / np.sqrt(2)
        return self

    def CNOT(self, control, target):
        """ apply CNOT gate with control and target qubits """
        old_state = copy.deepcopy(self.state) # avoid overwriting

        for i,_ in enumerate(old_state):
            if i & (1 << (control - 1)): # check if control qubit is 1
                
                index = i ^ (1 << (target - 1))
                self.state[i] = old_state[index]
                self.state[index] = old_state[i]
        return self




if __name__ == "__main__":
    print("Task 4:")
    x4 = quant(2)
    print("Initial state:")
    x4.print_state(precision=7)
    x4.H(2).X(1)
    print("After applying H_2 and X_1:")
    x4.print_state(precision=7)
    
    print("############ what was to show ############")
    print("Task 5:")
    x5 = quant(2)
    print("Initial state:")
    x5.print_state(precision=7)
    print("After applying H_1 and CNOT_(1,2):")
    x5.H(1).CNOT(1,2)
    x5.print_state(precision=7)
    print("############################################") 

    print("More complex example:")
    print("Initial state:")
    x = quant(3, coeffs=range(1,9))
    x.print_state(precision=7)
    
    x = x.H(1).CNOT(1,2).CNOT(2,3).H(3)
    print("After applying H_1, CNOT_(1,2), CNOT_(2,3), and H_3:")
    x.print_state(precision=7)



    
