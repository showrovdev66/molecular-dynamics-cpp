# Project 1: 2D Lennard-Jones Molecular Dynamics

## Why I made this project

I wanted to understand the basic structure of a molecular dynamics program before using larger simulation packages.

Instead of starting with a large codebase, I used only 36 particles in two dimensions. This made it easier to follow every part of the calculation.

## What the program does

The C++ program:

1. places 36 particles on a `6 x 6` grid;
2. gives them random starting velocities using a fixed random seed;
3. removes the average velocity so the complete system does not drift;
4. calculates Lennard-Jones forces between particle pairs;
5. applies periodic boundary conditions with the minimum-image convention;
6. moves the particles using the velocity-Verlet algorithm;
7. saves energy data to `results/energy.csv`;
8. saves particle coordinates to `results/trajectory.xyz`.

The Python script reads those two output files and creates `results/md_summary.png`.

## Main simulation settings

| Setting | Value |
|---|---:|
| Number of particles | 36 |
| Box size | 12 x 12 |
| Time step | 0.005 |
| Number of steps | 4000 |
| LJ epsilon | 1.0 |
| LJ sigma | 1.0 |
| Cutoff | 2.5 sigma |
| Saved frame interval | 20 steps |
| Random seed | 42 |

All quantities are in reduced Lennard-Jones units.

## How I did it

For two particles separated by distance `r`, the Lennard-Jones potential is

```text
U(r) = 4 epsilon [(sigma/r)^12 - (sigma/r)^6]
```

The short-range part is strongly repulsive and the longer-range part is attractive.

I calculate each pair only once. The force added to particle `i` is subtracted from particle `j`, which follows Newton's third law.

The simulation uses a square periodic box. If a particle leaves one side, it comes back through the opposite side. For force calculations I use the nearest periodic image of each particle.

The velocity-Verlet method updates position first, recalculates the forces, and then finishes the velocity update using the old and new forces.

## Example result

The included example run produced 200 saved frames.

![Energy and final particle positions](results/md_summary.png)

The example energy file starts with total energy about `8.007` and ends at about `7.794`. This is approximately a `2.7%` change over the saved run.

I therefore use the total-energy curve as a basic numerical check, but I do not describe the model as perfectly energy conserving. A finite time step and the simple hard cutoff can affect the energy behavior.

## Files

```text
01-lennard-jones-md/
├── README.md
├── src/
│   └── lj_md_simulation.cpp
├── scripts/
│   └── plot_trajectory.py
└── results/
    ├── energy.csv
    ├── trajectory.xyz
    └── md_summary.png
```

## Build and run

From this project folder:

```bash
g++ -O2 -std=c++17 -Wall -Wextra -pedantic src/lj_md_simulation.cpp -o lj_md_simulation
./lj_md_simulation
```

On Windows PowerShell with MinGW/MSYS2:

```powershell
g++ -O2 -std=c++17 -Wall -Wextra -pedantic src/lj_md_simulation.cpp -o lj_md_simulation.exe
.\lj_md_simulation.exe
```

Then create the plot:

```bash
python scripts/plot_trajectory.py
```

## Python packages

The plotting script uses:

- pandas
- matplotlib

Install them from the repository root with:

```bash
python -m pip install -r requirements.txt
```

## References

- J. E. Lennard-Jones, *Proc. R. Soc. A* **106**, 463-477 (1924). DOI: [10.1098/rspa.1924.0082](https://doi.org/10.1098/rspa.1924.0082)
- L. Verlet, *Physical Review* **159**, 98-103 (1967). DOI: [10.1103/PhysRev.159.98](https://doi.org/10.1103/PhysRev.159.98)
- W. C. Swope et al., *J. Chem. Phys.* **76**, 637-649 (1982). DOI: [10.1063/1.442716](https://doi.org/10.1063/1.442716)
- M. P. Allen and D. J. Tildesley, *Computer Simulation of Liquids*, 2nd ed., Oxford University Press (2017).
