// Cell-list optimization benchmark for a 2D Lennard-Jones force calculation
//
// The program compares two ways of finding interacting particle pairs:
// 1. a direct all-pairs search, and
// 2. a cell-list search that checks only nearby spatial cells.
//
// The cell-list method is approximately O(N) at fixed density and fixed cutoff
// because the average number of particles in nearby cells stays bounded as the
// system becomes larger.

#include <algorithm>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

struct Particle {
    double x = 0.0;
    double y = 0.0;
    double fx = 0.0;
    double fy = 0.0;
};

constexpr double EPSILON = 1.0;
constexpr double SIGMA = 1.0;
constexpr double CUTOFF = 2.5 * SIGMA;
constexpr double DENSITY = 0.5;

// Return the shortest displacement in a periodic box.
double minimum_image(double distance, double box_length) {
    if (distance > 0.5 * box_length) {
        distance -= box_length;
    }
    if (distance < -0.5 * box_length) {
        distance += box_length;
    }
    return distance;
}

// Put particles on a square grid at a fixed number density.
double initialize_grid(std::vector<Particle>& particles, int particles_per_side) {
    const int n = particles_per_side * particles_per_side;
    const double box_length = std::sqrt(n / DENSITY);
    const double spacing = box_length / particles_per_side;

    int index = 0;
    for (int i = 0; i < particles_per_side; ++i) {
        for (int j = 0; j < particles_per_side; ++j) {
            particles[index].x = (i + 0.5) * spacing;
            particles[index].y = (j + 0.5) * spacing;
            ++index;
        }
    }

    return box_length;
}

// Direct method: test every unique particle pair.
double compute_forces_naive(std::vector<Particle>& particles, double box_length) {
    for (auto& particle : particles) {
        particle.fx = 0.0;
        particle.fy = 0.0;
    }

    const double sigma6 = std::pow(SIGMA, 6);
    const double sigma12 = sigma6 * sigma6;
    const double cutoff2 = CUTOFF * CUTOFF;
    double potential_energy = 0.0;

    const int n = static_cast<int>(particles.size());
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            const double dx = minimum_image(particles[i].x - particles[j].x, box_length);
            const double dy = minimum_image(particles[i].y - particles[j].y, box_length);
            const double r2 = dx * dx + dy * dy;

            if (r2 >= cutoff2 || r2 <= 1.0e-12) {
                continue;
            }

            const double r2_inv = 1.0 / r2;
            const double r6_inv = r2_inv * r2_inv * r2_inv;
            const double force_factor =
                24.0 * EPSILON * r2_inv * r6_inv *
                (2.0 * sigma12 * r6_inv - sigma6);

            particles[i].fx += force_factor * dx;
            particles[i].fy += force_factor * dy;
            particles[j].fx -= force_factor * dx;
            particles[j].fy -= force_factor * dy;

            potential_energy +=
                4.0 * EPSILON *
                (sigma12 * r6_inv * r6_inv - sigma6 * r6_inv);
        }
    }

    return potential_energy;
}

// Cell-list method: place particles into spatial buckets and only inspect the
// current cell plus its eight neighboring cells.
double compute_forces_cell_list(std::vector<Particle>& particles, double box_length) {
    for (auto& particle : particles) {
        particle.fx = 0.0;
        particle.fy = 0.0;
    }

    const double sigma6 = std::pow(SIGMA, 6);
    const double sigma12 = sigma6 * sigma6;
    const double cutoff2 = CUTOFF * CUTOFF;
    double potential_energy = 0.0;

    const int n = static_cast<int>(particles.size());
    const int cells_per_side = std::max(3, static_cast<int>(box_length / CUTOFF));
    const double cell_size = box_length / cells_per_side;

    std::vector<std::vector<int>> cells(cells_per_side * cells_per_side);
    std::vector<int> particle_cell_x(n);
    std::vector<int> particle_cell_y(n);

    auto cell_coordinate = [cell_size, cells_per_side](double coordinate) {
        int cell = static_cast<int>(coordinate / cell_size);
        cell = std::max(0, cell);
        cell = std::min(cells_per_side - 1, cell);
        return cell;
    };

    // Build the cell list.
    for (int i = 0; i < n; ++i) {
        const int cx = cell_coordinate(particles[i].x);
        const int cy = cell_coordinate(particles[i].y);
        particle_cell_x[i] = cx;
        particle_cell_y[i] = cy;
        cells[cx * cells_per_side + cy].push_back(i);
    }

    // Search only nearby cells. Periodic wrapping is applied to cell indices.
    for (int i = 0; i < n; ++i) {
        const int cx = particle_cell_x[i];
        const int cy = particle_cell_y[i];

        for (int dx_cell = -1; dx_cell <= 1; ++dx_cell) {
            for (int dy_cell = -1; dy_cell <= 1; ++dy_cell) {
                const int neighbor_x =
                    ((cx + dx_cell) % cells_per_side + cells_per_side) % cells_per_side;
                const int neighbor_y =
                    ((cy + dy_cell) % cells_per_side + cells_per_side) % cells_per_side;

                for (int j : cells[neighbor_x * cells_per_side + neighbor_y]) {
                    if (j <= i) {
                        continue;  // Each pair is evaluated only once.
                    }

                    const double dx =
                        minimum_image(particles[i].x - particles[j].x, box_length);
                    const double dy =
                        minimum_image(particles[i].y - particles[j].y, box_length);
                    const double r2 = dx * dx + dy * dy;

                    if (r2 >= cutoff2 || r2 <= 1.0e-12) {
                        continue;
                    }

                    const double r2_inv = 1.0 / r2;
                    const double r6_inv = r2_inv * r2_inv * r2_inv;
                    const double force_factor =
                        24.0 * EPSILON * r2_inv * r6_inv *
                        (2.0 * sigma12 * r6_inv - sigma6);

                    particles[i].fx += force_factor * dx;
                    particles[i].fy += force_factor * dy;
                    particles[j].fx -= force_factor * dx;
                    particles[j].fy -= force_factor * dy;

                    potential_energy +=
                        4.0 * EPSILON *
                        (sigma12 * r6_inv * r6_inv - sigma6 * r6_inv);
                }
            }
        }
    }

    return potential_energy;
}

// Compare both methods on exactly the same particle positions before timing them.
bool verify_correctness() {
    constexpr int particles_per_side = 20;
    constexpr int n = particles_per_side * particles_per_side;

    std::vector<Particle> naive_particles(n);
    std::vector<Particle> cell_particles(n);

    const double box_length = initialize_grid(naive_particles, particles_per_side);
    initialize_grid(cell_particles, particles_per_side);

    const double naive_energy = compute_forces_naive(naive_particles, box_length);
    const double cell_energy = compute_forces_cell_list(cell_particles, box_length);

    double maximum_force_difference = 0.0;
    for (int i = 0; i < n; ++i) {
        maximum_force_difference = std::max(
            maximum_force_difference,
            std::abs(naive_particles[i].fx - cell_particles[i].fx));
        maximum_force_difference = std::max(
            maximum_force_difference,
            std::abs(naive_particles[i].fy - cell_particles[i].fy));
    }

    const double energy_difference = std::abs(naive_energy - cell_energy);
    const bool passed =
        maximum_force_difference < 1.0e-9 && energy_difference < 1.0e-9;

    std::cout << "=== Correctness check (N=" << n << ") ===\n";
    std::cout << "Maximum force component difference: "
              << maximum_force_difference << '\n';
    std::cout << "Potential energy difference: " << energy_difference << '\n';
    std::cout << "Result: " << (passed ? "PASS" : "FAIL") << "\n\n";

    return passed;
}

void run_benchmark() {
    const std::vector<int> particles_per_side = {10, 20, 30, 40, 50, 70};
    constexpr int repeats = 15;

    std::filesystem::create_directories("results");
    std::ofstream csv_file("results/benchmark_results.csv");
    if (!csv_file) {
        std::cerr << "Could not create results/benchmark_results.csv\n";
        return;
    }

    csv_file << "N,naive_ms,cell_list_ms,speedup\n";

    std::cout << "=== Benchmark: direct pairs vs. cell list ===\n";
    std::cout << std::left << std::setw(10) << "N"
              << std::setw(18) << "Direct (ms)"
              << std::setw(20) << "Cell list (ms)"
              << std::setw(10) << "Speedup" << '\n';

    for (int side : particles_per_side) {
        const int n = side * side;
        std::vector<Particle> particles(n);
        const double box_length = initialize_grid(particles, side);

        const auto direct_start = std::chrono::high_resolution_clock::now();
        for (int repeat = 0; repeat < repeats; ++repeat) {
            compute_forces_naive(particles, box_length);
        }
        const auto direct_end = std::chrono::high_resolution_clock::now();

        const auto cell_start = std::chrono::high_resolution_clock::now();
        for (int repeat = 0; repeat < repeats; ++repeat) {
            compute_forces_cell_list(particles, box_length);
        }
        const auto cell_end = std::chrono::high_resolution_clock::now();

        const double direct_ms =
            std::chrono::duration<double, std::milli>(direct_end - direct_start).count() /
            repeats;
        const double cell_ms =
            std::chrono::duration<double, std::milli>(cell_end - cell_start).count() /
            repeats;
        const double speedup = direct_ms / cell_ms;

        std::cout << std::fixed << std::setprecision(4)
                  << std::left << std::setw(10) << n
                  << std::setw(18) << direct_ms
                  << std::setw(20) << cell_ms
                  << std::setw(10) << speedup << '\n';

        csv_file << n << ',' << direct_ms << ',' << cell_ms << ',' << speedup << '\n';
    }

    std::cout << "\nSaved benchmark data to results/benchmark_results.csv\n";
}

int main() {
    if (!verify_correctness()) {
        std::cerr << "The two force methods do not agree, so the benchmark was stopped.\n";
        return 1;
    }

    run_benchmark();
    return 0;
}
