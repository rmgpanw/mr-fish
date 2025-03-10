# MR Fish: Proteomics and Mendelian Randomisation Analysis

Code to reproduce figures and tables from our publication: "Harnessing confounding and genetic pleiotropy to identify causes of disease through proteomics and Mendelian randomisation – 'MR Fish'"

## System Requirements

### Software Dependencies
- R (>= 4.1.0)
- RStudio (recommended for development)
- Docker (optional, for containerized execution)
- R package dependencies are detailed in `renv.lock`

### Operating System
- Tested on macOS (>= 10.15)
- Should work on Linux and Windows with appropriate R installation

## Installation Guide

1. Create a GitHub Codespace to run the analysis in your browser (5 minutes), or clone this repository locally:
```bash
git clone https://github.com/username/mr-fish.git
cd mr-fish
```

2. Install required R packages (5-10 minutes):
```R
# Install renv if not already installed
install.packages("renv")

# Initialize the project environment and install dependencies
renv::restore()
```

3. Run the analysis (2-3 minutes):
```R
targets::tar_make()
```

### Expected Outputs

The analysis will generate the following files in the `output/` directory:
- Figure 3.png
- Figure 4.png
- Figure S1.png
- Figure S2.png
- Figure S3.png
- Figure S4.png
- Figure A1.png
- sTable1.html
- sTable2.html
- sTable3.html
- sTable4.html

Expected runtime: ~5 minutes on a standard desktop computer

## License

This project is licensed under the MIT License. See the [LICENSE.md](LICENSE.md) file for details.
