# MobiDL 2

![plateforms:Linux-64](https://img.shields.io/badge/plateform-Linux--64-green?&style=for-the-badge)

MobiDL 2 is a collection of tools wrapped in WDL to be used in any WDL pipelines.

## __WARNING : BETA VERSION DO NOT USE IN PRODUCTION__

## Introduction

This repo provide a set of workflow in WDL (see [workflows](#workflows-implemented))
and a set of tools wrapped in WDL tasks (see [tools](tools.md)).

The code follows WDL specifications as much as possible ([WDL-spec](https://github.com/openwdl/wdl/blob/main/versions/1.0/SPEC.md)).

## Installation

This repo include an easy install by conda ([check the installation guide for conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html)).

1. Clone the repo:

`$ git clone https://github.com/mobidic/MobiDL2.0.git`

2. Install the dependencies:

`$ conda env create -f environment.yml`

3. Load the environnement:

`$ conda activate MobiDL2`

> If you have any issue with the installation please report the issue [here](https://github.com/mobidic/MobiDL2.0/issues/new?assignees=Char-Al&labels=%3Asnake%3A+Installation+bug&template=installation-issue-with-conda.md&title=Installation+with+conda+%3A+%5Bdetails%5D).

## Workflows

### Subworkflows

- postProcessAlignment
  - treatment of alignment files (sort, markdup, baseRecalibrator, compression)
- postProcessVCF
  - process a VCF to apply filter and normalization
