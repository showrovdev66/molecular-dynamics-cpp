# Project 2: Cell-List Optimization Benchmark

## Why I made this project

In the first molecular dynamics project, the force calculation checks every possible particle pair.

That method is simple and useful for learning, but it becomes expensive when the number of particles increases. I made this project to understand one common way MD programs reduce the number of unnecessary pair checks.

## The main idea

The direct method checks roughly

```text
N(N - 1) / 2
```

unique pairs, so its cost grows approximately as `O(N^2)`.

The optimized method divides the periodic box into square cells. The cell size is at least as large as the interaction cutoff.

A particle can then interact only with particles in:

- its own cell;
- the eight neighboring cells.

At fixed density and fixed cutoff, the average number of particles in those nearby cells stays roughly limited as the total system grows. For that reason the practical scaling approaches approximately `O(N)`.

## Why I included a correctness test

A faster program is not useful if it calculates the wrong forces.

Before measuring speed, the program creates the same 400-particle configuration for both methods and compares:

- every force component;
- total potential energy.

The benchmark only continues if the differences are below `1e-9`.

This was important for me because it gave me a simple rule: verify an optimization against a slower reference implementation before trusting benchmark numbers.

## Benchmark sizes

The program tests square systems with:

```text
N = 100, 400, 900, 1600, 2500, 4900
```

Each force calculation is repeated 15 times and the average time is reported.

The output is also saved to:

```text
results/benchmark_results.csv
```

Timing values depend strongly on CPU, compiler, optimization flags, and background processes. The important result is the scaling trend, not one exact speedup number.

## Files

```text
02-cell-list-benchmark/
├── README.md
├── src/
│   └── neighbor_list_benchmark.cpp
└── results/
    └── benchmark_results.csv   # created after running the program
```

## Build and run

Linux/macOS with GCC:

```bash
g++ -O2 -std=c++17 -Wall -Wextra -pedantic src/neighbor_list_benchmark.cpp -o neighbor_list_benchmark
./neighbor_list_benchmark
```

Windows PowerShell with MinGW/MSYS2:

```powershell
g++ -O2 -std=c++17 -Wall -Wextra -pedantic src/neighbor_list_benchmark.cpp -o neighbor_list_benchmark.exe
.\neighbor_list_benchmark.exe
```

## What I learned

This project helped me understand:

- algorithmic scaling;
- spatial decomposition;
- periodic cell indexing;
- why a cutoff makes local searching possible;
- why numerical optimization should be checked for correctness first;
- how benchmarking can separate algorithm quality from just code readability.

## References

- M. P. Allen and D. J. Tildesley, *Computer Simulation of Liquids*, 2nd ed., Oxford University Press (2017).
- J. E. Lennard-Jones, *Proc. R. Soc. A* **106**, 463-477 (1924). DOI: [10.1098/rspa.1924.0082](https://doi.org/10.1098/rspa.1924.0082)
