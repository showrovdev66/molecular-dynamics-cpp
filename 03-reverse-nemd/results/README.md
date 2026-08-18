# Reverse-NEMD Results

This folder contains the numerical output from my completed reverse-NEMD simulation.

## Files

### `temperature_profile.csv`

This file contains the final average temperature profile across the 20 spatial slabs used in the simulation.

The profile is used to measure the temperature gradient needed for the thermal-conductivity calculation.

### `nemd_summary.txt`

This file contains the main numerical results and validation checks from the simulation.

The completed run gave:

- measured thermal conductivity: `kappa* = 7.181061`
- published comparison value: `kappa* = 7.100000`
- relative difference: `1.141705%`
- momentum drift: `0.000000`
- energy drift: `0.001922%`
- left branch `R^2`: `0.989686`
- right branch `R^2`: `0.970907`

The small energy drift and high `R^2` values indicate that the short simulation produced a stable temperature profile suitable for this learning exercise.

The comparison value is used only as a reference. This project is not intended to claim an exact reproduction of the published simulation.