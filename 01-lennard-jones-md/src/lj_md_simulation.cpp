// 2D Lennard-Jones molecular dynamics simulator
//
// This program is a small learning project that shows the basic parts of a
// molecular dynamics (MD) simulation: particle initialization, force
// calculation, periodic boundaries, time integration, and result output.
//
// Main references:
// J. E. Lennard-Jones, Proc. R. Soc. A 106, 463-477 (1924).
// DOI: 10.1098/rspa.1924.0082
// W. C. Swope et al., J. Chem. Phys. 76, 637-649 (1982).
// DOI: 10.1063/1.442716

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <vector>

// Reduced Lennard-Jones units are used, so epsilon, sigma, and mass are 1.
constexpr double EPSILON = 1.0;
constexpr double SIGMA = 1.0;
constexpr double MASS = 1.0;
constexpr double BOX_LENGTH = 12.0;
constexpr double DT = 0.005;
constexpr double CUTOFF = 2.5 * SIGMA;

constexpr int N_PER_SIDE = 6;
constexpr int N_PARTICLES = N_PER_SIDE * N_PER_SIDE;
constexpr int N_STEPS = 4000;
constexpr int DUMP_EVERY = 20;

struct Particle {
    double x = 0.0;
    double y = 0.0;
    double vx = 0.0;
    double vy = 0.0;
    double fx = 0.0;
    double fy = 0.0;
};

// Return the shortest periodic displacement between two particles.
double minimum_image(double distance, double box_length) {
    if (distance > 0.5 * box_length) {
        distance -= box_length;
    }
    if (distance < -0.5 * box_length) {
        distance += box_length;
    }
    return distance;
}

// Put a coordinate back inside the periodic box [0, L).
double wrap_coordinate(double coordinate, double box_length) {
    while (coordinate < 0.0) {
        coordinate += box_length;
    }
    while (coordinate >= box_length) {
        coordinate -= box_length;
    }
    return coordinate;
}

// Calculate Lennard-Jones forces for every unique particle pair.
// The function also returns the total potential energy.
double compute_forces(std::vector<Particle>& particles) {
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
            const double dx = minimum_image(particles[i].x - particles[j].x, BOX_LENGTH);
            const double dy = minimum_image(particles[i].y - particles[j].y, BOX_LENGTH);
            const double r2 = dx * dx + dy * dy;

            // Ignore pairs outside the cutoff and protect against r = 0.
            if (r2 >= cutoff2 || r2 <= 1.0e-12) {
                continue;
            }

            const double r2_inv = 1.0 / r2;
            const double r6_inv = r2_inv * r2_inv * r2_inv;

            const double force_factor =
                24.0 * EPSILON * r2_inv * r6_inv *
                (2.0 * sigma12 * r6_inv - sigma6);

            // Newton's third law lets one pair calculation update both particles.
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

// Start particles on a regular grid and give them reproducible random velocities.
void initialize_particles(std::vector<Particle>& particles) {
    std::mt19937 random_engine(42);
    std::uniform_real_distribution<double> velocity_distribution(-1.0, 1.0);

    const double spacing = BOX_LENGTH / N_PER_SIDE;
    int index = 0;

    for (int i = 0; i < N_PER_SIDE; ++i) {
        for (int j = 0; j < N_PER_SIDE; ++j) {
            Particle particle;
            particle.x = (i + 0.5) * spacing;
            particle.y = (j + 0.5) * spacing;
            particle.vx = velocity_distribution(random_engine);
            particle.vy = velocity_distribution(random_engine);
            particles[index++] = particle;
        }
    }

    // Remove center-of-mass motion so the whole system does not drift.
    double vx_sum = 0.0;
    double vy_sum = 0.0;

    for (const auto& particle : particles) {
        vx_sum += particle.vx;
        vy_sum += particle.vy;
    }

    const double vx_mean = vx_sum / particles.size();
    const double vy_mean = vy_sum / particles.size();

    for (auto& particle : particles) {
        particle.vx -= vx_mean;
        particle.vy -= vy_mean;
    }
}

double kinetic_energy(const std::vector<Particle>& particles) {
    double energy = 0.0;

    for (const auto& particle : particles) {
        energy += 0.5 * MASS *
                  (particle.vx * particle.vx + particle.vy * particle.vy);
    }

    return energy;
}

int main() {
    std::vector<Particle> particles(N_PARTICLES);
    initialize_particles(particles);

    double potential_energy = compute_forces(particles);

    // Keep generated data separate from source code.
    std::filesystem::create_directories("results");
    std::ofstream trajectory_file("results/trajectory.xyz");
    std::ofstream energy_file("results/energy.csv");

    if (!trajectory_file || !energy_file) {
        std::cerr << "Could not create files inside the results folder.\n";
        return 1;
    }

    energy_file << "step,kinetic,potential,total\n";

    for (int step = 0; step < N_STEPS; ++step) {
        // Velocity-Verlet step 1: update positions using the current forces.
        for (auto& particle : particles) {
            particle.x = wrap_coordinate(
                particle.x + particle.vx * DT +
                    0.5 * (particle.fx / MASS) * DT * DT,
                BOX_LENGTH);

            particle.y = wrap_coordinate(
                particle.y + particle.vy * DT +
                    0.5 * (particle.fy / MASS) * DT * DT,
                BOX_LENGTH);
        }

        // Save the old force because velocity-Verlet uses old and new forces.
        std::vector<double> old_fx(N_PARTICLES);
        std::vector<double> old_fy(N_PARTICLES);

        for (int i = 0; i < N_PARTICLES; ++i) {
            old_fx[i] = particles[i].fx;
            old_fy[i] = particles[i].fy;
        }

        potential_energy = compute_forces(particles);

        // Velocity-Verlet step 2: finish the velocity update.
        for (int i = 0; i < N_PARTICLES; ++i) {
            particles[i].vx +=
                0.5 * (old_fx[i] + particles[i].fx) / MASS * DT;
            particles[i].vy +=
                0.5 * (old_fy[i] + particles[i].fy) / MASS * DT;
        }

        // Save one frame and one energy sample every DUMP_EVERY steps.
        if (step % DUMP_EVERY == 0) {
            const double ke = kinetic_energy(particles);
            const double total_energy = ke + potential_energy;

            energy_file << step << ',' << ke << ',' << potential_energy << ','
                        << total_energy << '\n';

            trajectory_file << N_PARTICLES << '\n';
            trajectory_file << "step " << step << '\n';

            for (const auto& particle : particles) {
                trajectory_file << "Ar " << particle.x << ' ' << particle.y
                                << " 0.0\n";
            }
        }

        if (step % 500 == 0) {
            const double ke = kinetic_energy(particles);
            std::cout << "step " << std::setw(5) << step
                      << "  KE = " << std::fixed << std::setprecision(4) << ke
                      << "  PE = " << potential_energy
                      << "  Total E = " << (ke + potential_energy) << '\n';
        }
    }

    std::cout << "\nFinished. Results were written to the results folder.\n";
    return 0;
}
