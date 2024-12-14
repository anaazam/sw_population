# Estimating the Population of Sex Workers

This repository contains data and replication code for estimating prostitution population dynamics using the single-registration capture-recapture approach. The methods are implemented in both Stata and R.

## Directory Structure

```
.
├── data/
│   ├── raw/       # Original Stata `.dta` files
│   ├── processed/ # Converted and cleaned `.csv` files
├── scripts/
│   ├── stata/     # Original Stata `.do` scripts
│   ├── r/         # Translated R scripts
├── results/       # Output files or tables from analysis
└── README.md      # This file
```

## Files

### Data
- `agg_dummies_13.dta`: Aggregated dummy variables for analysis.
- `belgium_dummies.dta`: Data specific to the Belgium model.
- `netherlands_dummies.dta`: Data specific to the Netherlands model.

### Scripts
#### Stata Scripts
- `pop_estimate.do`: Main estimation model.
- `belgium.do`: Country-specific model for Belgium.
- `netherlands.do`: Country-specific model for the Netherlands.

#### R Scripts
- `pop_estimate.R`: Translated main estimation model.
- `belgium_model.R`: Translated Belgium-specific model.
- `netherlands_model.R`: Translated Netherlands-specific model.

## Installation and Requirements

### R Dependencies
Install the required R packages before running the scripts:

```R
install.packages(c("haven", "dplyr", "tidyverse", "readr"))
```

### Conversion Notes
- Stata `.dta` files are converted to `.csv` for compatibility.
- Scripts are re-implemented in R for reproducibility.

## Usage

1. Clone the repository:
   ```bash
   git clone https://github.com/anaazam/sw_population.git
   cd sw_population
   ```

2. Process the data:
   - Raw Stata files are stored in `data/raw/`.
   - Converted `.csv` files are saved in `data/processed/`.

3. Run analyses:
   - Use R scripts in `scripts/r/` to replicate results.

## Methodology

The methodology is based on the paper "Estimating the prostitution population using online data: a single-registration approach for the Netherlands and Belgium" by Azam et al. (2024). The main model employs the Zelterman estimator to derive population counts using capture-recapture techniques.

## License
This project is licensed under the MIT License.
