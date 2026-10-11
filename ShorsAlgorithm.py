from qiskit import QuantumCircuit
import numpy as np
import math
import random
import matplotlib.pyplot as plt

last_qc = None

N = int(input("Input desired composite: "))
factors = set()

def generate_primes(limit):
    primes = []
    for num in range(2, limit + 1):
        is_prime = True
        for i in range(2, int(math.sqrt(num)) + 1):
            if num % i == 0:
                is_prime = False
                break
        if is_prime:
            primes.append(num)
    return primes

primes_to_test = generate_primes(N // 2)

def Shors(a):
    binary = bin(N)[2:][::-1]
    L = len(binary)
    print(f"Simulating an {L} qubit Shor's Algorithm for integer {N}!")

    qc = QuantumCircuit(4*L+2, 2*L)

    def QFT(x):
        start_idx = x * L
        end_idx = start_idx + L - 1
        for i in range(end_idx-start_idx+1):
            qc.h(start_idx+i)
            theta = np.pi/2
            for j in range(end_idx-start_idx-i):
                qc.cp(theta, start_idx+i+1+j, start_idx+i)
                theta = theta/2
        for i in range(int((end_idx-start_idx+1)/2)):
            qc.swap(start_idx+i, end_idx-i)

    def IQFT(x):
        start_idx = x * L
        end_idx = start_idx + L - 1
        for i in range(int((end_idx-start_idx+1)/2)):
            qc.swap(start_idx+i, end_idx-i)
        for i in range(end_idx-start_idx, -1, -1):
            for j in range(end_idx-start_idx-i-1, -1, -1):
                theta = -(np.pi/(2**(j+1)))
                qc.cp(theta, start_idx+i+1+j, start_idx+i)
            qc.h(start_idx+i)

    def IQFT_2L():
        for i in range(int(2*L/2)):
            qc.swap(i, 2*L-1-i)
        for i in range(2*L-1, -1, -1):
            for j in range(2*L-1-i-1, -1, -1):
                theta = -(np.pi/(2**(j+1)))
                qc.cp(theta, i+1+j, i)
            qc.h(i)

    def DraperAdder(x):
        QFT(x)
        start_idx = x * L
        end_idx = start_idx + L - 1
        for i in range(end_idx-start_idx+1):
            for j in range(2*L):
                theta = (np.pi/(2**(j-i)))
                qc.cp(theta, j, start_idx+i)
        IQFT(x)

    def InverseDraperAdder(x):
        QFT(x)
        start_idx = x * L
        end_idx = start_idx + L - 1
        for i in range(end_idx-start_idx+1):
            for j in range(2*L):
                theta = -(np.pi/(2**(j-i)))
                qc.cp(theta, j, start_idx+i)
        IQFT(x)

    def ModularAdder(x, control, val_string):
        start_idx = x * L
        end_idx = start_idx + L - 1
        QFT(x)
        for i in range(end_idx-start_idx+1):
            for j, char in enumerate(val_string):
                if int(char) == 1:
                    theta = (np.pi/(2**(j-i)))
                    qc.cp(theta, control, start_idx+i)
        for i in range(end_idx-start_idx+1):
            for j, char in enumerate(binary):
                if int(char) == 1:
                    theta = -(np.pi/(2**(j-i)))
                    qc.p(theta, start_idx+i)
        IQFT(x)
        qc.cx(end_idx, 4*L)
        QFT(x)
        for i in range(end_idx-start_idx+1):
            for j, char in enumerate(binary):
                if int(char) == 1:
                    theta = (np.pi/(2**(j-i)))
                    qc.cp(theta, 4*L, start_idx+i)
        IQFT(x)
        QFT(x)
        for i in range(end_idx-start_idx+1):
            for j, char in enumerate(val_string):
                if int(char) == 1:
                    theta = -(np.pi/(2**(j-i)))
                    qc.cp(theta, control, start_idx+i)
        IQFT(x)
        qc.x(end_idx)
        qc.cx(end_idx, 4*L)
        qc.x(end_idx)
        QFT(x)
        for i in range(end_idx-start_idx+1):
            for j, char in enumerate(val_string):
                if int(char) == 1:
                    theta = (np.pi/(2**(j-i)))
                    qc.cp(theta, control, start_idx+i)
        IQFT(x)

    def ModularMultiplier():
        for i in range(2*L):
            val = (a * (2**i)) % N
            val_bin = bin(val)[2:].zfill(L)[::-1]
            ModularAdder(3, i, val_bin)
            
        for i in range(L):
            qc.swap(2*L+i, 3*L+i)
            
        a_inv = pow(a, -1, N)
        for i in range(2*L):
            val = (a_inv * (2**i)) % N
            val_bin = bin(val)[2:].zfill(L)[::-1]
            start_idx = 3 * L
            end_idx = start_idx + L - 1
            QFT(3)
            for k in range(end_idx-start_idx+1):
                for j, char in enumerate(val_bin):
                    if int(char) == 1:
                        theta = -(np.pi/(2**(j-k)))
                        qc.cp(theta, i, start_idx+k)
            IQFT(3)

    # 1. Apply hadamards to every qubit in top register
    for i in range(2*L):
        qc.h(i)

    # 2. Prepare lower register with 1 (Little-Endian)
    qc.x(2*L)

    qc.barrier()

    # 3. Calculate a^x = f(x) mod(N)
        # a. Draper Adder
    DraperAdder(2)

    qc.barrier()

        # b. Modular Adder (to calculate mod(N))

        # c. Modular Multiplier
    ModularMultiplier()
    IQFT_2L()
    qc.measure(range(2*L), range(2*L))
    global last_qc
    last_qc = qc

for a in primes_to_test:
    print(f"\n--- Testing base a = {a} ---")
    g = math.gcd(a, N)
    
    if g > 1:
        print(f"Found a factor classically via GCD: {g}")
        factors.add(g)
        factors.add(N // g)
        # 2. The 'break' is removed so it keeps testing the next 'a'!
        
    else:
        print(f"a and N are coprime. Building quantum circuit for a={a}...")
        Shors(a)

# 3. Print the final collected set of factors
print(f"\nAll discovered factors for {N}: {sorted(list(factors))}")

# if N <= 50:
    # fig = last_qc.draw(output='mpl', fold=60)
    # fig.set_size_inches(24, 12)
    # fig.suptitle(f"Shor's Algorithm Circuit (N={N}, a={a})")
    # plt.show()