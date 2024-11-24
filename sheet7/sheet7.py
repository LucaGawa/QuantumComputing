import numpy as np

def GCD(x,y):
    """
        calculate the greatest common divisor (GCD) of the two integers x and y using the Euclidean algorithm.
    """
    if x == 0:
        return y
    else:
        return GCD(y % x, x)

def extendedGCD(a,b):
    if b == 0:
        return a, 1, 0
    else:
        d, s, t = extendedGCD(b, a % b)
        return d, t, s - (a // b) * t

def orderFinder(x, N):
    """
        find least positive integer r that fullfills
        x^r = 1(mod N) returns r    
    """
    r = 1
    while (x**r % N) != 1:
        r += 1
    else:
        return r


def shore(N):
    """
        factorize the integer N by using Shore's algorithm.
    """
    i=0
    while True:
        x = np.random.randint(1, N)
        x = GCD(x, N)
        i += 1
        if x != 1:
            return x
        else:
            r = orderFinder(x,N)
            if not ((r%2 == 1) or (x**(r/2)%N == -1)): # r odd or x^r/2 =-1 (mod n)
                return GCD(x**(r/2)-1,N)

def genM(string):
    """
        generate the integer M from a given string
    """
    M = 0
    l = len(string)
    for i,c in enumerate(string):
        # map a -> 0, b -> 1, c -> 2, ...
        M += (ord(c)-97)*26**(l-i-1) # -97 since a is 97 in ascii
    return M

def genString(M):
    """
        generate the string from the integer M
    """
    l =  4
    string = ""
    for i in range(l):
        string += chr((M // (26**(l-i-1)))%26 + 97)

    return string

def modularInverse(e, phi):
    """
        calculate the modular inverse of e mod phi
    """
    gcd, s, t = extendedGCD(e, phi)
    if gcd != 1:
        raise ValueError("modular inverse does not exist")
    return (s % phi + phi) % phi

def mod_pow(base, exponent, modulus):
    result = 1
    base = base % modulus
    while exponent > 0:
        if exponent % 2 == 1:
            result = (result * base) % modulus
        exponent = exponent >> 1
        base = (base * base) % modulus
    return result

def decrypt(N, e, string):
    q = shore(N)
    p = int(N/q)
    phi = (p-1)*(q-1)
    _,s,_ = extendedGCD(e,phi)
    d = modularInverse(e,phi)
    M = genM(string)
    Md = mod_pow(M,d,N)
    return genString(Md)

if __name__ == "__main__":
    N = 71219839 
    e = 4187899
    string = "bwtlqv"
    print(decrypt(N,e, string))

    N = 121130231
    e = 20579
    string = "blhhay"
    print(decrypt(N,e, string))





