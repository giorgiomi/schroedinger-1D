# Schrödinger equation in 1D

The equation we want to solve is the following:

```math
    i\frac{\partial}{\partial t}\psi(x,t) = \left[-\frac{1}{2}\frac{\partial^2}{\partial x^2} + V(x)\right]\psi(x,t)\, ,
```

where $\hbar = m = 1$. We are interested in the region $x \in [-L,L]$ and we take as boundary conditions $\psi(-L,t) = \psi(L,t) = 0$ for all $t$.

## Compilation and execution
### Single simulation
In order to run a single simulation, compile in the repository directory using:
```
    make eu
    make cn
    make ts
```
depending on the solver you want to use (`eu` stands for Euler, `cn` for Crank-Nicolson and `ts` for Trotter-Suzuki). An executable named `<method>.x` will be generated. Use
```
    ./<method>.x <N> <dt> <V0> <a>
```
to execute it, where `N` is the number of spatial divisions, `dt` is the time step, `V0` is the potential strength parameter, `a` is the potential well position parameter. Then, to plot the results one can use
```
    python3 plots/probability.py
    python3 plots/motion.py
    python3 plots/energy.py
```

### Multiple simulation
In order to study the behaviour of the system with different parameters, one can run multiple simulations. This is automated via some `.sh` scripts, named
```
    multi_a.sh
    multi_V.sh
    multi_N.sh
```
Their usage is explained when launching `./multi_<a/V/N>.sh`. The compilation is included in these scripts, so one does not need to run `make` before. Some data visualization is also provided.