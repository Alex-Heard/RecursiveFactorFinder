# Shor's Algorithm & Quantum Arithmetic Circuit Generator

An end-to-end implementation of **Shor's Algorithm** built in Python with **Qiskit**. This program decomposes composite integers by testing potential coprimes and synthesizing full gate-level quantum arithmetic circuits for modular exponentiation and order finding.

---

## Technical Overview

Rather than utilizing abstract unitaries or black-box oracles, this implementation builds quantum arithmetic directly from foundational gate primitives:

- **Quantum Fourier Transform (QFT / IQFT):** Implemented via parameterized controlled-phase rotations ($CP$) and bit-reversal SWAP networks.
- **Draper Adder:** Performs out-of-place addition in the Fourier frequency basis using phase kickback without carry-propagation overhead.
- **Modular Adder:** Constructs controlled addition modulo $N$ using Draper stages, bit-flag auxiliary registers, and uncomputation.
- **Modular Multiplier:** Computes controlled modular multiplication ($a \cdot x \pmod N$) and reverses ancillary register entanglement using modular inverse arithmetic ($a^{-1} \pmod N$)[cite: 6].
- **Classical Co-processing:** Generates prime candidates below $N/2$, checks for classical gcd factors, and constructs quantum circuits whenever $a$ and $N$ are coprime[cite: 6].

---

## Circuit Architecture

For an integer $N$ with bit-length $L = \lceil\log_2 N\rceil$[cite: 6]:
- **Register Allocation:** $4L + 2$ qubits total[cite: 6]
  - $2L$ Counting / Control qubits (initialized into uniform superposition via Hadamard gates)[cite: 6]
  - Work and target registers for modular arithmetic subcircuits[cite: 6]
  - Auxiliary carry/flag qubits for modular comparison[cite: 6]
- **Measurement:** $2L$ classical bits measuring the phase register following an inverse QFT ($IQFT_{2L}$)[cite: 6].

---

## Prerequisites

```bash
pip install qiskit numpy matplotlib
