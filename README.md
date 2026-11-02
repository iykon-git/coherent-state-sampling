# Coherent-state sampling

Python implementation of the coherent-state sampling (CSS) algorithm described in **“Coherent state sampling for molecular vibronic spectra”**.

Associated preprint: https://doi.org/10.21203/rs.3.rs-8593810/v1

Project repository: https://github.com/iykon-git/coherent-state-sampling

## Description

The script samples vibrational configurations within the linear coupling model using vibrational wavenumbers and Huang–Rhys factors as input.

For each mode $i$, the vibrational occupation number is independently sampled from a Poisson distribution with mean $S_i$, where $S_i$ is the Huang–Rhys factor:

$$
n_i \sim \operatorname{Poisson}(S_i).
$$

One CSS event is a complete configuration $\mathbf{n}=(n_1,\ldots,n_M)$ of all $M$ input modes. Sampling is performed without an occupation-number cutoff.

For each observed vibrational state, the script provides its energy, relative intensity estimated from its sampled probability, and line assignment specified by the number of vibrational quanta in each mode.

## Requirements

- Python 3
- NumPy

The script was tested with **Python 3.11.14** and **NumPy 2.4.2**.

Install NumPy, if needed:

```bash
python -m pip install numpy
```

## Files

- `coherent_state_sampling.py` — CSS implementation.
- `anthracene_input.txt` — example input containing the anthracene data used in the associated article.
- `README.md` — documentation.
- `LICENSE` — MIT License.

## Input data

The input is a whitespace-separated text file with exactly two columns and one row per vibrational mode:

| Column | Quantity | Unit |
| --- | --- | --- |
| 1 | Vibrational wavenumber | cm⁻¹ |
| 2 | Huang–Rhys factor | Dimensionless |

Example format:

```text
# wavenumber_cm^-1  Huang_Rhys_factor
500.0              0.10
1000.0             0.20
1500.0             0.05
```

The numbers above illustrate the format only.

Lines beginning with `#` are comments. Both columns must contain finite values, and Huang–Rhys factors must be non-negative. The input row order defines the mode numbering in the output.

The supplied `anthracene_input.txt` contains the vibrational wavenumbers and Huang–Rhys factors used for anthracene in the associated article. All modes are retained, including those with very small Huang–Rhys factors. Users may select a subset of modes by editing the input file.

## Configuration

Place the input file in the same directory as the script. Edit the constants in the `USER SETTINGS` section before running:

```python
MOLECULE_NAME = "anthracene"
INPUT_FILE_NAME = "anthracene_input.txt"

N_EVENTS = 100_000
BLOCK_SIZE = 100_000
RNG_SEED = 101

ROUND_ENERGIES = True
OUTPUT_DIR_NAME = None
```

| Setting | Description |
| --- | --- |
| `MOLECULE_NAME` | Label used in the run summary and default output folder name. |
| `INPUT_FILE_NAME` | Input filename, resolved relative to the script directory. |
| `N_EVENTS` | Total number of full multimode CSS events. Must be positive. |
| `BLOCK_SIZE` | Maximum number of events generated per block. Must be positive. |
| `RNG_SEED` | Integer random seed; `None` requests a fresh random initialization. |
| `ROUND_ENERGIES` | If `True`, round input wavenumbers to the nearest integer using `numpy.rint` before calculating state energies. |
| `OUTPUT_DIR_NAME` | Output folder name, relative to the script directory. If `None`, a name is generated from the molecule label and seed. |

The Huang–Rhys factors determine the sampled occupation numbers. The vibrational wavenumbers determine the corresponding state energies.

## Usage

From the directory containing the script, run:

```bash
python coherent_state_sampling.py
```

Alternatively, open the script in Spyder and run it after editing the settings.

Progress is printed after each sampling block.

## Output

With the settings above, the script creates the directory:

```text
css_basic_anthracene_seed_101/
```

It contains two files.

### `css_unique_states.csv`

Each row represents a unique observed vibrational configuration.

| Column | Description |
| --- | --- |
| `rank` | Row number, starting at 1. |
| `n_mode_1` … `n_mode_M` | Number of vibrational quanta in each mode, in input row order; these specify the line assignment. |
| `energy_cm^-1` | Vibrational excitation energy expressed in cm⁻¹. |
| `count` | Number of times the configuration was sampled. |
| `sampled_probability` | Relative frequency of the configuration, providing an estimate of its normalized line intensity. |

For a configuration $\mathbf{n}$, the reported vibrational excitation energy in wavenumber units is

$$
\tilde{\nu}_{\mathbf{n}}=\sum_{i=1}^{M}n_i\tilde{\nu}_i,
$$

where $\tilde{\nu}_i$ is the input wavenumber after optional rounding. The zero-occupation configuration defines the energy origin.

The sampled probability is

$$
\hat{p}_{\mathbf{n}}=\frac{C_{\mathbf{n}}}{R},
$$

where $C_{\mathbf{n}}$ is the configuration count and $R$ is `N_EVENTS`.

Rows are sorted by decreasing count, with equal counts sorted by increasing energy. Distinct configurations remain separate even when their energies coincide. Only observed configurations appear in the table.

### `run_summary.json`

Records:

- Molecule label and input filename.
- Number of vibrational modes.
- Number of CSS events.
- Block size and random seed.
- Wavenumber rounding setting.
- Number of unique observed states.
- Sum of configuration counts.

Existing output files in the selected output directory are overwritten when the script is run again.

## Reproducibility

To reproduce a run, retain the input file and its row order, the script version, all settings, and the Python and NumPy versions.

Use an integer `RNG_SEED` and the same `BLOCK_SIZE` and software environment. Samples are generated mode by mode within each block, so changing the block size can change the sampled results even with the same seed.

Block processing limits the size of the temporary event array. Counts for all unique observed configurations are retained in memory throughout the run.

## Author

Iurii Konyshev

## Citation

If you use this code in your research, please cite the associated preprint:

**Coherent state sampling for molecular vibronic spectra**  
https://doi.org/10.21203/rs.3.rs-8593810/v1

Please also cite the Zenodo archive of the software release used in your work. Release archives can be located through the GitHub repository:

https://github.com/iykon-git/coherent-state-sampling

## License

This code is distributed under the MIT License. See the `LICENSE` file for details.
