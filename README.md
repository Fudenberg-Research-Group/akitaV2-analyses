## Description

This repository houses scripts and notebooks used to generate and analyze data with the AkitaV2 model.

Scripts used for cross-species AkitaV2 training and model weights are available from the [Basenji repository](https://github.com/calico/basenji/tree/master/manuscripts/akita/v2).

Preprint available here: *link to be added very soon*

## Prerequisites

This repository depends on `akita_utils` available here: [akita_utils repository](https://github.com/Fudenberg-Research-Group/akita_utils).
Please install `akita_utils` before running notebooks or tutorials in this repository. 

## Contents:

### 1. experiments/
Contains all experiments conducted during the project. Refer to individual experiment folders for detailed documentation and results. Each experiment is essentially a unique (sequence generator, TSV table) pair. The TSV file specifies the genomic positions and experiment parameters, while the sequence generator specifies how to construct the sequence for input into AkitaV2.

### 2. figures/
Includes default plotting configurations and example maps utilized in illustrations.

### 3. tutorials/
Includes two notebooks exlaining and visualizing disruption and insertion experiments.

### 4. input data/
Consists of data preprocessing scripts and crucial TSV tables containing information about CTCF sites, including their genomic locations. The CTCF sites analyzed do not overlap with each other or with repeatable elements.

## Contact Information

Feedback and questions are appreciated. Please contact us at: fudenber at usc fullstop edu & smaruj at usc fullstop edu.

We welcome contributions from the community. Please fork this repository and submit pull requests. For major changes, please open an issue first to discuss what you would like to change.

## License

This project is licensed under the MIT License - see the LICENSE file for details.
