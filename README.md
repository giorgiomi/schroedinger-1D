# Schrödinger equation in 1D

The equation we want to solve is the following:

```math
    i\frac{\partial}{\partial t}\psi(x,t) = \left[-\frac{1}{2}\frac{\partial^2}{\partial x^2} + V(x)\right]\psi(x,t)\, ,
```

where $\hbar = m = 1$. We are interested in the region $x \in [-L,L]$ and we take as boundary conditions $\psi(-L,t) = \psi(L,t) = 0$ for all $t$.

## Compilation and execution
Compile in the repository directory using:
```
    make cn
    make eu
```
depending on the solver you want to use (`cn` stands for Crank-Nicolson and `eu` for Euler). An executable named `run` will be generated. Use
```
    ./run
```
to execute it. Then, to plot the results one can use
```
    python3 plots/normalization.py
    python3 plots/motion.py
```