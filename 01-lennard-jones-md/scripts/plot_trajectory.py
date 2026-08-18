"""Plot the output from the 2D Lennard-Jones C++ simulation.

The script reads the energy history and XYZ trajectory from the project's
results folder. It creates one summary image containing an energy plot and the
final particle positions.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


BOX_LENGTH = 12.0
PROJECT_DIR = Path(__file__).resolve().parent.parent
RESULTS_DIR = PROJECT_DIR / "results"
ENERGY_FILE = RESULTS_DIR / "energy.csv"
TRAJECTORY_FILE = RESULTS_DIR / "trajectory.xyz"
OUTPUT_IMAGE = RESULTS_DIR / "md_summary.png"


def read_xyz(filename: Path) -> list[list[tuple[float, float]]]:
    """Read a simple XYZ trajectory and return its 2D particle positions."""
    frames = []

    with filename.open("r", encoding="utf-8") as file:
        while True:
            particle_count_line = file.readline()
            if not particle_count_line:
                break

            particle_count = int(particle_count_line.strip())
            file.readline()  # The second line stores the simulation step.

            frame = []
            for _ in range(particle_count):
                parts = file.readline().split()
                x = float(parts[1])
                y = float(parts[2])
                frame.append((x, y))

            frames.append(frame)

    return frames


def main() -> None:
    if not ENERGY_FILE.exists() or not TRAJECTORY_FILE.exists():
        raise FileNotFoundError(
            "Run the C++ simulation first so energy.csv and trajectory.xyz exist."
        )

    energy = pd.read_csv(ENERGY_FILE)
    frames = read_xyz(TRAJECTORY_FILE)

    if not frames:
        raise ValueError("The trajectory file does not contain any frames.")

    figure, axes = plt.subplots(1, 2, figsize=(12, 5))

    # Left side: check how kinetic, potential, and total energy change.
    axes[0].plot(energy["step"], energy["kinetic"], label="Kinetic")
    axes[0].plot(energy["step"], energy["potential"], label="Potential")
    axes[0].plot(energy["step"], energy["total"], label="Total", linewidth=2)
    axes[0].set_xlabel("Step")
    axes[0].set_ylabel("Energy")
    axes[0].set_title("Energy during the simulation")
    axes[0].legend()

    # Right side: show where the particles are in the final saved frame.
    final_frame = frames[-1]
    x_positions = [position[0] for position in final_frame]
    y_positions = [position[1] for position in final_frame]

    axes[1].scatter(x_positions, y_positions, s=80, edgecolors="black")
    axes[1].set_xlim(0, BOX_LENGTH)
    axes[1].set_ylim(0, BOX_LENGTH)
    axes[1].set_aspect("equal")
    axes[1].set_xlabel("x")
    axes[1].set_ylabel("y")
    axes[1].set_title(f"Final saved frame ({len(frames) - 1})")

    figure.tight_layout()
    figure.savefig(OUTPUT_IMAGE, dpi=150)

    print(f"Read {len(frames)} trajectory frames.")
    print(f"Saved: {OUTPUT_IMAGE}")

    plt.show()


if __name__ == "__main__":
    main()
