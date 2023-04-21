Calculate PSI from the sequence
================

Program to  calculate Protein Stability Index from the sequence. ML models were developed based on experimental stability datasets for a 24/23-mer covering the N-/C-terminus of the human proteome using the CatBoost regressor. The performance of the final models was evaluated using the testing set and an R2 coefficient, reaching the values of 0.796/0.812 for the N-terminus with initiator methionine cleaved/not cleaved, respectively, and 0.815 for the C-terminus (the highest possible value of R2 coefficient is 1). See the paper for details and the [degronopedia tutorial](https://degronopedia.com/degronopedia/tutorial#ML) for more information.

The web version of this tool (and more!) is available at: [https://degronopedia.com/](https://degronopedia.com/).

[![action status](https://github.com/filipsPL/degronopedia-ml-psi/actions/workflows/thefirst.yml/badge.svg)](https://github.com/filipsPL/degronopedia-ml-psi/actions/workflows/thefirst.yml)

- [Calculate PSI from the sequence](#calculate-psi-from-the-sequence)
  - [Set up](#set-up)
  - [Usage](#usage)
  - [Interpretation](#interpretation)
  - [Authors](#authors)
  - [How to cite](#how-to-cite)

## Set up

```sh
# clone the repo
git clone --depth 1 https://github.com/filipsPL/degronopedia-ml-psi

# create a new conda environment
conda env create -f conda.yml`
```

## Usage

Input is a sequence in a plain text format, eg:

```text
MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGPDEAPRMPEAAPPVAPAPAAPTPAAPAPAPSWPLSSSVPSQKTYQGSYGFRLGFLHSGTAKSVTCTYSPALNKMFCQLAKTCPVQLWVDSTPPPGTRVRAMAIYKQSQHMTEVVRRCPHHERCSDSDGLAPPQHLIRVEGNLRVEYLDDRNTFRHSVVVPYEPPEVGSDCTTIHYNYMCNSSCMGGMNRRPILTIITLEDSSGNLLGRNSFEVRVCACPGRDRRT
```

Options:

- the file with sequence
- which terminus calculate PSI for. Choices:
  - `C` for C-terminus
  - `NiMetNo` for N-terminus with initiator Met cleaved
  - `NiMetYes` for N-terminus with initiator Met **NOT** cleaved

Running the program:

```sh
# activate the environment
conda activate dp

# run the program
./calculate-desc.py --sequence sequence.txt --type NiMetYes
```

Output will be:

```text
N-terminus with initiator Met NOT cleaved
[2023-04-21 15:21:00] Predicting...
[2023-04-21 15:21:00] Predicting done.
Predicted PSI: 5.24
```

## Interpretation

Predictions are made based on the datasets of experimental PSI values, which describe the stability of protein N-/C-terminus in an artificial system where 23-mers covering the termini of nearly entire human proteome were conjugated to GFP protein, and their stability was measured relative to the stability of DsRed protein translated from the same transcript using the Global Protein Stability (GPS) high-throughput technique (Koren et al., 2018 and Timms et al., 2019). Therefore, these values provide insight into the stability of the N-/C-terminus of the query but to a limited extent. Several peptides with low PSI values were experimentally validated to be degraded by the cullin-RING E3 ligase complexes by the authors of the aforementioned GPS studies. However, medium or higher PSI values do not rule out the regulation of such termini by N-/C-degron pathways, as other factors may influence this, including tissue specificity, posttranslational modifications, stress conditions, etc.

## Authors

Szulc, N. A., Stefaniak, F.

## How to cite

Szulc, N. A., Stefaniak, F., Piechota, M., Cappannini, A., Bujnicki, J. M., & Pokrzywa, W. (2022). DEGRONOPEDIA - a web server for proteome-wide inspection of degrons. bioRxiv. <https://doi.org/10.1101/2022.05.19.492622>
