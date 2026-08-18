// 3D reverse-NEMD thermal conductivity calculation for a Lennard-Jones fluid
//
// This is an educational implementation of the Müller-Plathe reverse
// non-equilibrium molecular dynamics method. The simulation imposes a heat
// flux by exchanging velocities between two slabs, measures the temperature
// gradient, and estimates thermal conductivity from Fourier's law.
//
// Main reference:
// F. Müller-Plathe, J. Chem. Phys. 106, 6082-6085 (1997).
// DOI: 10.1063/1.473271

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <string>
#include <vector>

struct Particle {
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
    double vx = 0.0;
    double vy = 0.0;
    double vz = 0.0;
    double fx = 0.0;
    double fy = 0.0;
    double fz = 0.0;
};

// Lennard-Jones reduced units.
constexpr double EPSILON = 1.0;
constexpr double SIGMA = 1.0;
constexpr double MASS = 1.0;
constexpr double CUTOFF = 3.0 * SIGMA;
constexpr double RHO_STAR = 0.849;
constexpr double T_STAR = 0.7;
constexpr double DT = 0.005;

// This smaller educational system is not an exact copy of the published
// 2592-particle validation cell. It uses the same reduced density, temperature,
// and cutoff as a useful reference state point.
constexpr int N_SIDE = 12;
constexpr int N_PARTICLES = N_SIDE * N_SIDE * N_SIDE;

constexpr int N_EQUIL = 3000;
constexpr int RESCALE_EVERY = 50;
constexpr int N_PRODUCTION = 12000;
constexpr int N_SWAP = 20;
constexpr int N_SLABS = 20;
constexpr int SAMPLE_EVERY = 5;
constexpr int DASHBOARD_EVERY = 2000;

constexpr double REFERENCE_KAPPA = 7.1;

struct LinearFit {
    double slope = 0.0;
    double r_squared = 0.0;
};

double minimum_image(double distance, double box_length) {
    if (distance > 0.5 * box_length) {
        distance -= box_length;
    }
    if (distance < -0.5 * box_length) {
        distance += box_length;
    }
    return distance;
}

double wrap_coordinate(double coordinate, double box_length) {
    while (coordinate < 0.0) {
        coordinate += box_length;
    }
    while (coordinate >= box_length) {
        coordinate -= box_length;
    }
    return coordinate;
}

int slab_index(double z, double box_length) {
    int slab = static_cast<int>(z / box_length * N_SLABS);
    slab = std::max(0, slab);
    slab = std::min(N_SLABS - 1, slab);
    return slab;
}

// Calculate LJ forces with a 3D cell list. Each particle only searches its own
// cell and the 26 neighboring cells.
double compute_forces(std::vector<Particle>& particles, double box_length) {
    for (auto& particle : particles) {
        particle.fx = 0.0;
        particle.fy = 0.0;
        particle.fz = 0.0;
    }

    const double sigma6 = std::pow(SIGMA, 6);
    const double sigma12 = sigma6 * sigma6;
    const double cutoff2 = CUTOFF * CUTOFF;
    double potential_energy = 0.0;

    const int n = static_cast<int>(particles.size());
    const int cells_per_side = std::max(3, static_cast<int>(box_length / CUTOFF));
    const double cell_size = box_length / cells_per_side;

    static std::vector<std::vector<int>> cells;
    const std::size_t number_of_cells =
        static_cast<std::size_t>(cells_per_side * cells_per_side * cells_per_side);

    if (cells.size() != number_of_cells) {
        cells.assign(number_of_cells, {});
    } else {
        for (auto& cell : cells) {
            cell.clear();
        }
    }

    auto get_cell = [cell_size, cells_per_side](double coordinate) {
        int cell = static_cast<int>(coordinate / cell_size);
        cell = std::max(0, cell);
        cell = std::min(cells_per_side - 1, cell);
        return cell;
    };

    std::vector<int> cell_x(n);
    std::vector<int> cell_y(n);
    std::vector<int> cell_z(n);

    for (int i = 0; i < n; ++i) {
        const int cx = get_cell(particles[i].x);
        const int cy = get_cell(particles[i].y);
        const int cz = get_cell(particles[i].z);

        cell_x[i] = cx;
        cell_y[i] = cy;
        cell_z[i] = cz;
        cells[(cx * cells_per_side + cy) * cells_per_side + cz].push_back(i);
    }

    for (int i = 0; i < n; ++i) {
        for (int dx_cell = -1; dx_cell <= 1; ++dx_cell) {
            for (int dy_cell = -1; dy_cell <= 1; ++dy_cell) {
                for (int dz_cell = -1; dz_cell <= 1; ++dz_cell) {
                    const int nx =
                        ((cell_x[i] + dx_cell) % cells_per_side + cells_per_side) %
                        cells_per_side;
                    const int ny =
                        ((cell_y[i] + dy_cell) % cells_per_side + cells_per_side) %
                        cells_per_side;
                    const int nz =
                        ((cell_z[i] + dz_cell) % cells_per_side + cells_per_side) %
                        cells_per_side;

                    const auto& neighbor_cell =
                        cells[(nx * cells_per_side + ny) * cells_per_side + nz];

                    for (int j : neighbor_cell) {
                        if (j <= i) {
                            continue;
                        }

                        const double rx =
                            minimum_image(particles[i].x - particles[j].x, box_length);
                        const double ry =
                            minimum_image(particles[i].y - particles[j].y, box_length);
                        const double rz =
                            minimum_image(particles[i].z - particles[j].z, box_length);
                        const double r2 = rx * rx + ry * ry + rz * rz;

                        if (r2 >= cutoff2 || r2 <= 1.0e-12) {
                            continue;
                        }

                        const double r2_inv = 1.0 / r2;
                        const double r6_inv = r2_inv * r2_inv * r2_inv;
                        const double force_factor =
                            24.0 * EPSILON * r2_inv * r6_inv *
                            (2.0 * sigma12 * r6_inv - sigma6);

                        particles[i].fx += force_factor * rx;
                        particles[i].fy += force_factor * ry;
                        particles[i].fz += force_factor * rz;
                        particles[j].fx -= force_factor * rx;
                        particles[j].fy -= force_factor * ry;
                        particles[j].fz -= force_factor * rz;

                        potential_energy +=
                            4.0 * EPSILON *
                            (sigma12 * r6_inv * r6_inv - sigma6 * r6_inv);
                    }
                }
            }
        }
    }

    return potential_energy;
}

// Build the starting configuration and set the initial reduced temperature.
double initialize_particles(std::vector<Particle>& particles) {
    const double box_length = std::cbrt(N_PARTICLES / RHO_STAR);
    const double spacing = box_length / N_SIDE;

    std::mt19937 random_engine(7);
    std::normal_distribution<double> velocity_distribution(0.0, 1.0);

    int index = 0;
    double vx_sum = 0.0;
    double vy_sum = 0.0;
    double vz_sum = 0.0;

    for (int i = 0; i < N_SIDE; ++i) {
        for (int j = 0; j < N_SIDE; ++j) {
            for (int k = 0; k < N_SIDE; ++k) {
                Particle particle;
                particle.x = (i + 0.5) * spacing;
                particle.y = (j + 0.5) * spacing;
                particle.z = (k + 0.5) * spacing;
                particle.vx = velocity_distribution(random_engine);
                particle.vy = velocity_distribution(random_engine);
                particle.vz = velocity_distribution(random_engine);

                vx_sum += particle.vx;
                vy_sum += particle.vy;
                vz_sum += particle.vz;
                particles[index++] = particle;
            }
        }
    }

    const double vx_mean = vx_sum / N_PARTICLES;
    const double vy_mean = vy_sum / N_PARTICLES;
    const double vz_mean = vz_sum / N_PARTICLES;

    double kinetic_energy = 0.0;
    for (auto& particle : particles) {
        particle.vx -= vx_mean;
        particle.vy -= vy_mean;
        particle.vz -= vz_mean;

        kinetic_energy += 0.5 * MASS *
                          (particle.vx * particle.vx +
                           particle.vy * particle.vy +
                           particle.vz * particle.vz);
    }

    const double current_temperature =
        2.0 * kinetic_energy / (3.0 * N_PARTICLES);
    const double scale = std::sqrt(T_STAR / current_temperature);

    for (auto& particle : particles) {
        particle.vx *= scale;
        particle.vy *= scale;
        particle.vz *= scale;
    }

    return box_length;
}

double total_kinetic_energy(const std::vector<Particle>& particles) {
    double energy = 0.0;
    for (const auto& particle : particles) {
        energy += 0.5 * MASS *
                  (particle.vx * particle.vx +
                   particle.vy * particle.vy +
                   particle.vz * particle.vz);
    }
    return energy;
}

double system_temperature(const std::vector<Particle>& particles) {
    return 2.0 * total_kinetic_energy(particles) /
           (3.0 * particles.size());
}

void rescale_temperature(std::vector<Particle>& particles, double target_temperature) {
    const double current_temperature = system_temperature(particles);
    const double scale = std::sqrt(target_temperature / current_temperature);

    for (auto& particle : particles) {
        particle.vx *= scale;
        particle.vy *= scale;
        particle.vz *= scale;
    }
}

// Perform one velocity-Verlet integration step.
double velocity_verlet_step(std::vector<Particle>& particles, double box_length) {
    for (auto& particle : particles) {
        particle.x = wrap_coordinate(
            particle.x + particle.vx * DT + 0.5 * particle.fx / MASS * DT * DT,
            box_length);
        particle.y = wrap_coordinate(
            particle.y + particle.vy * DT + 0.5 * particle.fy / MASS * DT * DT,
            box_length);
        particle.z = wrap_coordinate(
            particle.z + particle.vz * DT + 0.5 * particle.fz / MASS * DT * DT,
            box_length);
    }

    std::vector<double> old_fx(particles.size());
    std::vector<double> old_fy(particles.size());
    std::vector<double> old_fz(particles.size());

    for (std::size_t i = 0; i < particles.size(); ++i) {
        old_fx[i] = particles[i].fx;
        old_fy[i] = particles[i].fy;
        old_fz[i] = particles[i].fz;
    }

    const double potential_energy = compute_forces(particles, box_length);

    for (std::size_t i = 0; i < particles.size(); ++i) {
        particles[i].vx += 0.5 * (old_fx[i] + particles[i].fx) / MASS * DT;
        particles[i].vy += 0.5 * (old_fy[i] + particles[i].fy) / MASS * DT;
        particles[i].vz += 0.5 * (old_fz[i] + particles[i].fz) / MASS * DT;
    }

    return potential_energy;
}

// Exchange the hottest particle in the cold slab with the coldest particle in
// the hot slab. Equal masses mean the velocity-vector swap conserves total
// kinetic energy and total momentum.
double muller_plathe_swap(std::vector<Particle>& particles, double box_length) {
    int hottest_in_cold = -1;
    int coldest_in_hot = -1;
    double largest_cold_ke = -1.0;
    double smallest_hot_ke = 1.0e300;

    for (std::size_t i = 0; i < particles.size(); ++i) {
        const int slab = slab_index(particles[i].z, box_length);
        const double kinetic_energy =
            0.5 * MASS *
            (particles[i].vx * particles[i].vx +
             particles[i].vy * particles[i].vy +
             particles[i].vz * particles[i].vz);

        if (slab == 0 && kinetic_energy > largest_cold_ke) {
            largest_cold_ke = kinetic_energy;
            hottest_in_cold = static_cast<int>(i);
        }

        if (slab == N_SLABS / 2 && kinetic_energy < smallest_hot_ke) {
            smallest_hot_ke = kinetic_energy;
            coldest_in_hot = static_cast<int>(i);
        }
    }

    if (hottest_in_cold < 0 || coldest_in_hot < 0) {
        return 0.0;
    }

    const double transferred_energy = largest_cold_ke - smallest_hot_ke;

    std::swap(particles[hottest_in_cold].vx, particles[coldest_in_hot].vx);
    std::swap(particles[hottest_in_cold].vy, particles[coldest_in_hot].vy);
    std::swap(particles[hottest_in_cold].vz, particles[coldest_in_hot].vz);

    return transferred_energy;
}

void print_temperature_profile(const std::vector<double>& temperatures, int step) {
    double minimum_temperature =
        *std::min_element(temperatures.begin(), temperatures.end());
    double maximum_temperature =
        *std::max_element(temperatures.begin(), temperatures.end());

    if (maximum_temperature - minimum_temperature < 1.0e-9) {
        maximum_temperature = minimum_temperature + 1.0e-9;
    }

    std::cout << "\nTemperature profile at production step " << step << '\n';

    for (int slab = 0; slab < N_SLABS; ++slab) {
        const int bar_length = static_cast<int>(
            40.0 * (temperatures[slab] - minimum_temperature) /
            (maximum_temperature - minimum_temperature));

        std::string label = "       ";
        if (slab == 0) {
            label = " [COLD]";
        } else if (slab == N_SLABS / 2) {
            label = " [HOT] ";
        }

        std::cout << "slab " << std::setw(2) << slab << label
                  << " T=" << std::fixed << std::setprecision(3)
                  << std::setw(6) << temperatures[slab]
                  << " |" << std::string(bar_length, '#') << '\n';
    }
}

LinearFit fit_line(const std::vector<double>& x, const std::vector<double>& y) {
    const int n = static_cast<int>(x.size());
    double sum_x = 0.0;
    double sum_y = 0.0;
    double sum_xx = 0.0;
    double sum_xy = 0.0;

    for (int i = 0; i < n; ++i) {
        sum_x += x[i];
        sum_y += y[i];
        sum_xx += x[i] * x[i];
        sum_xy += x[i] * y[i];
    }

    const double slope =
        (n * sum_xy - sum_x * sum_y) /
        (n * sum_xx - sum_x * sum_x);
    const double intercept = (sum_y - slope * sum_x) / n;
    const double mean_y = sum_y / n;

    double residual_sum = 0.0;
    double total_sum = 0.0;

    for (int i = 0; i < n; ++i) {
        const double predicted = slope * x[i] + intercept;
        residual_sum += (y[i] - predicted) * (y[i] - predicted);
        total_sum += (y[i] - mean_y) * (y[i] - mean_y);
    }

    return {slope, 1.0 - residual_sum / total_sum};
}

void write_results(
    const std::vector<double>& temperatures,
    double box_length,
    double kappa,
    double momentum_drift,
    double energy_drift_percent,
    const LinearFit& left_fit,
    const LinearFit& right_fit) {

    std::filesystem::create_directories("results");

    std::ofstream profile_file("results/temperature_profile.csv");
    profile_file << "slab,z_center,reduced_temperature\n";

    const double slab_thickness = box_length / N_SLABS;
    for (int slab = 0; slab < N_SLABS; ++slab) {
        const double z_center = (slab + 0.5) * slab_thickness;
        profile_file << slab << ',' << z_center << ',' << temperatures[slab] << '\n';
    }

    std::ofstream summary_file("results/nemd_summary.txt");
    summary_file << std::fixed << std::setprecision(6);
    summary_file << "Measured kappa*: " << kappa << '\n';
    summary_file << "Reference kappa*: " << REFERENCE_KAPPA << '\n';
    summary_file << "Relative difference (%): "
                 << std::abs(kappa - REFERENCE_KAPPA) / REFERENCE_KAPPA * 100.0 << '\n';
    summary_file << "Momentum drift: " << momentum_drift << '\n';
    summary_file << "Energy drift (%): " << energy_drift_percent << '\n';
    summary_file << "Left branch R^2: " << left_fit.r_squared << '\n';
    summary_file << "Right branch R^2: " << right_fit.r_squared << '\n';
}

int main() {
    std::vector<Particle> particles(N_PARTICLES);
    const double box_length = initialize_particles(particles);

    std::cout << "=== 3D reverse NEMD thermal conductivity ===\n";
    std::cout << "N=" << N_PARTICLES
              << "  L=" << box_length
              << "  rho*=" << RHO_STAR
              << "  T*=" << T_STAR
              << "  cutoff=" << CUTOFF << "\n\n";

    compute_forces(particles, box_length);

    // Equilibrate while periodically rescaling the temperature.
    std::cout << "Equilibrating...\n";
    for (int step = 0; step < N_EQUIL; ++step) {
        velocity_verlet_step(particles, box_length);

        if (step % RESCALE_EVERY == 0) {
            rescale_temperature(particles, T_STAR);
        }

        if (step % 1000 == 0) {
            std::cout << "equil step " << step
                      << "  T=" << system_temperature(particles) << '\n';
        }
    }

    double initial_px = 0.0;
    double initial_py = 0.0;
    double initial_pz = 0.0;
    for (const auto& particle : particles) {
        initial_px += particle.vx;
        initial_py += particle.vy;
        initial_pz += particle.vz;
    }

    const double initial_total_energy =
        total_kinetic_energy(particles) + compute_forces(particles, box_length);

    // Production run: no thermostat, only the Müller-Plathe velocity swaps.
    std::cout << "\nRunning production for " << N_PRODUCTION << " steps...\n";

    std::vector<double> temperature_sum(N_SLABS, 0.0);
    std::vector<int> sample_count(N_SLABS, 0);
    double total_transferred_energy = 0.0;

    for (int step = 0; step < N_PRODUCTION; ++step) {
        velocity_verlet_step(particles, box_length);

        if (step % N_SWAP == 0) {
            total_transferred_energy += muller_plathe_swap(particles, box_length);
        }

        if (step % SAMPLE_EVERY == 0) {
            std::vector<double> slab_kinetic_energy(N_SLABS, 0.0);
            std::vector<int> particles_in_slab(N_SLABS, 0);

            for (const auto& particle : particles) {
                const int slab = slab_index(particle.z, box_length);
                slab_kinetic_energy[slab] +=
                    0.5 * MASS *
                    (particle.vx * particle.vx +
                     particle.vy * particle.vy +
                     particle.vz * particle.vz);
                ++particles_in_slab[slab];
            }

            for (int slab = 0; slab < N_SLABS; ++slab) {
                if (particles_in_slab[slab] > 0) {
                    const double slab_temperature =
                        2.0 * slab_kinetic_energy[slab] /
                        (3.0 * particles_in_slab[slab]);
                    temperature_sum[slab] += slab_temperature;
                    ++sample_count[slab];
                }
            }
        }

        if (step > 0 && step % DASHBOARD_EVERY == 0) {
            std::vector<double> current_average(N_SLABS, 0.0);
            for (int slab = 0; slab < N_SLABS; ++slab) {
                current_average[slab] =
                    temperature_sum[slab] / std::max(1, sample_count[slab]);
            }
            print_temperature_profile(current_average, step);
        }
    }

    std::vector<double> average_temperature(N_SLABS, 0.0);
    for (int slab = 0; slab < N_SLABS; ++slab) {
        average_temperature[slab] =
            temperature_sum[slab] / std::max(1, sample_count[slab]);
    }
    print_temperature_profile(average_temperature, N_PRODUCTION);

    // Fit the two nearly linear branches, excluding the source/sink slabs.
    std::vector<double> left_z;
    std::vector<double> left_temperature;
    std::vector<double> right_z;
    std::vector<double> right_temperature;

    for (int slab = 1; slab < N_SLABS / 2; ++slab) {
        left_z.push_back(slab);
        left_temperature.push_back(average_temperature[slab]);
    }

    for (int slab = N_SLABS / 2 + 1; slab < N_SLABS; ++slab) {
        right_z.push_back(slab);
        right_temperature.push_back(average_temperature[slab]);
    }

    const LinearFit left_fit = fit_line(left_z, left_temperature);
    const LinearFit right_fit = fit_line(right_z, right_temperature);

    const double slab_thickness = box_length / N_SLABS;
    const double left_gradient = left_fit.slope / slab_thickness;
    const double right_gradient = -right_fit.slope / slab_thickness;
    const double average_gradient =
        0.5 * (std::abs(left_gradient) + std::abs(right_gradient));

    const double production_time = N_PRODUCTION * DT;
    const double area = box_length * box_length;

    // Periodic boundaries create two equal heat-flow paths, giving the factor 2.
    const double heat_flux =
        total_transferred_energy / (2.0 * area * production_time);
    const double kappa = heat_flux / average_gradient;

    double final_px = 0.0;
    double final_py = 0.0;
    double final_pz = 0.0;
    for (const auto& particle : particles) {
        final_px += particle.vx;
        final_py += particle.vy;
        final_pz += particle.vz;
    }

    const double momentum_drift = std::sqrt(
        (final_px - initial_px) * (final_px - initial_px) +
        (final_py - initial_py) * (final_py - initial_py) +
        (final_pz - initial_pz) * (final_pz - initial_pz));

    const double final_total_energy =
        total_kinetic_energy(particles) + compute_forces(particles, box_length);
    const double energy_drift_percent =
        std::abs(final_total_energy - initial_total_energy) /
        std::abs(initial_total_energy) * 100.0;

    std::cout << "\n=== Validation ===\n";
    std::cout << "Momentum drift: " << momentum_drift << '\n';
    std::cout << "Energy drift during production: " << energy_drift_percent << " %\n";
    std::cout << "Linear fit R^2: left=" << left_fit.r_squared
              << ", right=" << right_fit.r_squared << '\n';

    std::cout << "\n=== Result ===\n";
    std::cout << "Measured thermal conductivity kappa* = " << kappa << '\n';
    std::cout << "Published reference value kappa* = " << REFERENCE_KAPPA << '\n';
    std::cout << "Relative difference = "
              << std::abs(kappa - REFERENCE_KAPPA) / REFERENCE_KAPPA * 100.0
              << " %\n";

    write_results(
        average_temperature,
        box_length,
        kappa,
        momentum_drift,
        energy_drift_percent,
        left_fit,
        right_fit);

    std::cout << "Saved temperature_profile.csv and nemd_summary.txt in results/.\n";
    return 0;
}
