# -*- coding: utf-8 -*-
"""
Basic state-resolved Coherent-State Sampling (CSS).

Input: a two-column text file with vibrational energies and Huang-Rhys factors.

Output:
    1) css_unique_states.csv: observed states, energies, counts, sampled probabilities.
    2) run_summary.json: run settings and event-count summary.

States are sorted by decreasing sampled count.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Optional

import numpy as np


# ========================== USER SETTINGS ==========================

MOLECULE_NAME = "anthracene"

# Two columns: energy, Huang-Rhys factor.
INPUT_FILE_NAME = "anthracene_input.txt"

# Number of CSS events.
N_EVENTS = 100_000

# Events are generated in blocks to limit memory use.
BLOCK_SIZE = 100_000

# Integer for reproducible results; None for a different run each time.
RNG_SEED: Optional[int] = 101

# Round input vibrational energies to the nearest integer.
ROUND_ENERGIES = True

# Output folder. None gives css_basic_<molecule>_seed_<seed>.
OUTPUT_DIR_NAME: Optional[str] = None


# =========================== PATH SETUP ============================

try:
    BASE_DIR = Path(__file__).resolve().parent
except NameError:
    BASE_DIR = Path.cwd().resolve()

INPUT_FILE = BASE_DIR / INPUT_FILE_NAME

if OUTPUT_DIR_NAME is None:
    seed_label = "random" if RNG_SEED is None else str(RNG_SEED)
    OUTPUT_DIR_NAME = f"css_basic_{MOLECULE_NAME}_seed_{seed_label}"

OUT_DIR = BASE_DIR / OUTPUT_DIR_NAME



# =========================== INPUT DATA ============================

def load_modes(path: Path, round_energies: bool = True):
    """Load vibrational energies and Huang-Rhys factors."""
    data = np.loadtxt(path, dtype=float, comments="#")
    data = np.atleast_2d(data)

    if data.ndim != 2 or data.shape[1] != 2:
        raise ValueError(
            "The input file must contain exactly two columns: "
            "energy and Huang-Rhys factor. "
            f"Received shape {data.shape}."
        )

    energies = data[:, 0].astype(np.float64)
    hr_factors = data[:, 1].astype(np.float64)

    if round_energies:
        energies = np.rint(energies)

    if not np.all(np.isfinite(energies)):
        raise ValueError("All vibrational energies must be finite.")
    if not np.all(np.isfinite(hr_factors)):
        raise ValueError("All Huang-Rhys factors must be finite.")
    if np.any(hr_factors < 0):
        raise ValueError("Huang-Rhys factors must be non-negative.")

    return energies, hr_factors


# =========================== CSS SAMPLING ==========================

def generate_css_block(
    rng: np.random.Generator,
    hr_factors: np.ndarray,
    block_size: int,
):
    """
    Generate one unrestricted CSS block.

    For every mode i:
        n_i ~ Poisson(S_i)

    There is no occupation-number cutoff and no clamping.
    """
    n_modes = hr_factors.size
    block = np.empty((block_size, n_modes), dtype=np.int32)

    for mode_index, hr_factor in enumerate(hr_factors):
        samples = rng.poisson(hr_factor, size=block_size)

        if samples.size and samples.max() > np.iinfo(np.int32).max:
            raise OverflowError("A sampled occupation exceeds int32 range.")

        block[:, mode_index] = samples.astype(np.int32, copy=False)

    return block


def run_css(
    hr_factors: np.ndarray,
    n_events: int,
    block_size: int,
    seed: Optional[int],
):
    """Run CSS and return all unique observed states with their counts."""
    if n_events <= 0:
        raise ValueError("N_EVENTS must be positive.")
    if block_size <= 0:
        raise ValueError("BLOCK_SIZE must be positive.")

    rng = np.random.default_rng(seed)
    state_counts: dict[tuple[int, ...], int] = {}

    processed_events = 0

    while processed_events < n_events:
        current_block_size = min(block_size, n_events - processed_events)

        block = generate_css_block(
            rng=rng,
            hr_factors=hr_factors,
            block_size=current_block_size,
        )

        unique_states, counts = np.unique(
            block,
            axis=0,
            return_counts=True,
        )

        for state_array, count in zip(unique_states, counts):
            state_tuple = tuple(int(value) for value in state_array)
            state_counts[state_tuple] = (
                state_counts.get(state_tuple, 0) + int(count)
            )

        processed_events += current_block_size

        print(
            f"Processed {processed_events:>12,} / {n_events:,} events | "
            f"unique states: {len(state_counts):,}"
        )

    states = np.asarray(list(state_counts.keys()), dtype=np.int32)
    counts = np.asarray(list(state_counts.values()), dtype=np.int64)

    if int(counts.sum()) != n_events:
        raise RuntimeError(
            "Count conservation failed: "
            f"sum(counts)={int(counts.sum())}, N_EVENTS={n_events}."
        )

    return states, counts


# ============================= OUTPUT ==============================

def save_state_table(
    path: Path,
    states: np.ndarray,
    energies: np.ndarray,
    counts: np.ndarray,
    sampled_probabilities: np.ndarray,
):
    """Save observed states, energies, counts, and sampled probabilities."""
    total_energies = states.astype(np.float64) @ energies

    # Sort by decreasing count; for equal counts, by increasing energy.
    order = np.lexsort((total_energies, -counts))

    states = states[order]
    total_energies = total_energies[order]
    counts = counts[order]
    sampled_probabilities = sampled_probabilities[order]

    n_modes = states.shape[1]
    header = (
        ["rank"]
        + [f"n_mode_{i + 1}" for i in range(n_modes)]
        + ["energy_cm^-1", "count", "sampled_probability"]
    )

    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(header)

        for rank, row_index in enumerate(range(states.shape[0]), start=1):
            writer.writerow(
                [rank]
                + states[row_index].astype(int).tolist()
                + [
                    f"{total_energies[row_index]:.12g}",
                    int(counts[row_index]),
                    f"{sampled_probabilities[row_index]:.12e}",
                ]
            )


# =============================== MAIN ==============================

def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    print("Basic state-resolved Coherent-State Sampling")
    print(f"Input file: {INPUT_FILE}")
    print(f"N_EVENTS: {N_EVENTS:,} | BLOCK_SIZE: {BLOCK_SIZE:,} | RNG_SEED: {RNG_SEED}")

    energies, hr_factors = load_modes(
        INPUT_FILE,
        round_energies=ROUND_ENERGIES,
    )

    states, counts = run_css(
        hr_factors=hr_factors,
        n_events=N_EVENTS,
        block_size=BLOCK_SIZE,
        seed=RNG_SEED,
    )

    sampled_probabilities = counts.astype(np.float64) / float(N_EVENTS)

    output_table = OUT_DIR / "css_unique_states.csv"
    save_state_table(
        path=output_table,
        states=states,
        energies=energies,
        counts=counts,
        sampled_probabilities=sampled_probabilities,
    )

    summary = {
        "molecule_name": MOLECULE_NAME,
        "input_file": INPUT_FILE_NAME,
        "number_of_modes": int(hr_factors.size),
        "n_events": int(N_EVENTS),
        "block_size": int(BLOCK_SIZE),
        "rng_seed": RNG_SEED,
        "round_energies": bool(ROUND_ENERGIES),
        "number_of_unique_states": int(states.shape[0]),
        "sum_of_counts": int(counts.sum()),
    }

    summary_path = OUT_DIR / "run_summary.json"
    with open(summary_path, "w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2)

    print(f"Unique observed states: {states.shape[0]:,}")
    print(f"Saved state table: {output_table}")
    print(f"Saved summary: {summary_path}")
    print("Done.")


if __name__ == "__main__":
    main()
