# References

These are the main references behind the physical models and numerical methods used in this repository.

## Lennard-Jones interaction

1. J. E. Lennard-Jones, "On the Determination of Molecular Fields. II. From the Equation of State of a Gas," *Proceedings of the Royal Society A*, **106**, 463-477 (1924). DOI: [10.1098/rspa.1924.0082](https://doi.org/10.1098/rspa.1924.0082)

## Molecular-dynamics integration

2. L. Verlet, "Computer 'Experiments' on Classical Fluids. I. Thermodynamical Properties of Lennard-Jones Molecules," *Physical Review*, **159**, 98-103 (1967). DOI: [10.1103/PhysRev.159.98](https://doi.org/10.1103/PhysRev.159.98)

3. W. C. Swope, H. C. Andersen, P. H. Berens, and K. R. Wilson, "A computer simulation method for the calculation of equilibrium constants for the formation of physical clusters of molecules: Application to small water clusters," *The Journal of Chemical Physics*, **76**, 637-649 (1982). DOI: [10.1063/1.442716](https://doi.org/10.1063/1.442716)

## General molecular-simulation methods and cell lists

4. M. P. Allen and D. J. Tildesley, *Computer Simulation of Liquids*, 2nd ed., Oxford University Press (2017), ISBN 978-0-19-880319-5.

## Reverse non-equilibrium molecular dynamics

5. F. Müller-Plathe, "A simple nonequilibrium molecular dynamics method for calculating the thermal conductivity," *The Journal of Chemical Physics*, **106**, 6082-6085 (1997). DOI: [10.1063/1.473271](https://doi.org/10.1063/1.473271)

6. P. B. Allen and Y. Li, "Phonon thermal conductivity by non-local non-equilibrium molecular dynamics," arXiv:1412.3099. The Lennard-Jones-liquid test in this work reports `rho = 0.849/sigma^3`, `T = 0.7 epsilon/kB`, cutoff `3.0 sigma`, and `kappa = 7.1` in Lennard-Jones units when reproducing the Müller-Plathe calculation. [arXiv:1412.3099](https://arxiv.org/abs/1412.3099)

## Citation note

The source code in this repository is a learning implementation. The papers and book above are cited for the physical models and numerical methods. They are not claimed as sources of copied program code.
