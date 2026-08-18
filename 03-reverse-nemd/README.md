# Project 3: Reverse-NEMD Thermal Conductivity

## Why I made this project

After learning the basic MD loop and a cell-list optimization, I wanted to connect the same ideas to a transport-property calculation.

This project is a small 3D Lennard-Jones simulation based on the reverse non-equilibrium molecular dynamics method introduced by Florian Müller-Plathe.

## What reverse NEMD means

Normally, heat flows from a hot region to a cold region.

In the Müller-Plathe method, the simulation deliberately transfers kinetic energy in the opposite direction. The system responds by developing a temperature gradient.

The basic steps are:

1. divide the periodic simulation box into slabs;
2. choose one slab as the cold region and another slab half a box away as the hot region;
3. at fixed intervals, find the hottest particle in the cold slab;
4. find the coldest particle in the hot slab;
5. exchange their complete velocity vectors;
6. keep track of the kinetic energy moved by the exchanges;
7. measure the average temperature in each slab;
8. fit the two temperature-gradient branches;
9. calculate thermal conductivity using Fourier's law.

For equal-mass particles, exchanging complete velocity vectors preserves total kinetic energy and total momentum during the swap itself.

## Simulation settings

| Setting | Value |
|---|---:|
| Particles | 1728 (`12 x 12 x 12`) |
| Reduced density `rho*` | 0.849 |
| Reduced temperature `T*` | 0.7 |
| LJ cutoff | 3.0 sigma |
| Time step | 0.005 |
| Equilibration steps | 3000 |
| Production steps | 12000 |
| Number of slabs | 20 |
| Swap interval | 20 steps |
| Sampling interval | 5 steps |

## Important comparison note

The program uses a published Lennard-Jones state point as a reference, but this learning implementation is not an exact reproduction of the original published simulation.

The published comparison reported in later work used 2592 atoms, `rho = 0.849/sigma^3`, `T = 0.7 epsilon/kB`, periodic boundaries, and a `3.0 sigma` cutoff, and reported `kappa = 7.1` in Lennard-Jones units when reproducing the Müller-Plathe calculation.

My program uses fewer particles and a different simulation-cell geometry, so `7.1` should be treated as a reference value rather than an exact answer that this short run must reproduce.

## What the program checks

At the end it reports:

- measured thermal conductivity `kappa*`;
- difference from the reference value;
- momentum drift;
- total-energy drift during production;
- `R^2` for the two linear temperature-profile fits.

It also saves:

```text
results/temperature_profile.csv
results/nemd_summary.txt
```

## Build and run

Linux/macOS with GCC:

```bash
g++ -O2 -std=c++17 -Wall -Wextra -pedantic src/nemd_thermal_conductivity_3d.cpp -o nemd_thermal_conductivity_3d
./nemd_thermal_conductivity_3d
```

Windows PowerShell with MinGW/MSYS2:

```powershell
g++ -O2 -std=c++17 -Wall -Wextra -pedantic src/nemd_thermal_conductivity_3d.cpp -o nemd_thermal_conductivity_3d.exe
.\nemd_thermal_conductivity_3d.exe
```

This project is much more computationally expensive than the first two, so runtime depends on the computer.

## What I learned

This project helped me connect several ideas:

- 3D periodic molecular dynamics;
- extending a 2D cell list to 3D;
- equilibration and production stages;
- spatial temperature profiles;
- imposed heat flux;
- Fourier's law;
- linear regression and `R^2` as simple quality checks;
- the importance of system size and simulation length in transport calculations.

## References

- F. Müller-Plathe, "A simple nonequilibrium molecular dynamics method for calculating the thermal conductivity," *J. Chem. Phys.* **106**, 6082-6085 (1997). DOI: [10.1063/1.473271](https://doi.org/10.1063/1.473271)
- P. B. Allen and Y. Li, "Phonon thermal conductivity by non-local non-equilibrium molecular dynamics," arXiv:1412.3099. The Lennard-Jones test reports the reference state point and `kappa = 7.1` used here for comparison. [arXiv:1412.3099](https://arxiv.org/abs/1412.3099)
- M. P. Allen and D. J. Tildesley, *Computer Simulation of Liquids*, 2nd ed., Oxford University Press (2017).
