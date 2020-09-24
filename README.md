# MobiDL 2

![plateforms:Linux-64|osx-64](https://img.shields.io/badge/plateform-Linux--64|osx--64-green?&style=for-the-badge)

MobiDL 2 is a collection of tools wrapped in WDL to be used in any WDL pipelines.

## __WARNING : BETA VERSION DO NOT USE IN PRODUCTION__

## Introduction

This repo provide a set of workflow in WDL (see [workflows](#workflows-implemented))
and a set of tools wrapped in WDL tasks (see [tools](#tools-implemented)).

The code follows WDL specifications as much as possible ([WDL-spec](https://github.com/openwdl/wdl/blob/main/versions/1.0/SPEC.md)).

## Installation

This repo include an easy install by conda ([check the installation guide for conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html)).

1. Clone the repo:

`$ git clone https://github.com/mobidic/MobiDL2.0.git`

2. Install the dependencies:

`$ conda env create environment.yml`

3. Load the environnement:

`$ conda activate MobiDL2`

> If you have any issue with the installation please report the issue [here](https://github.com/mobidic/MobiDL2.0/issues/new?assignees=Char-Al&labels=%3Asnake%3A+Installation+bug&template=installation-issue-with-conda.md&title=Installation+with+conda+%3A+%5Bdetails%5D).

## Workflows implemented

### General

#### PrepareGenome

This workflow is developped to prepare a fasta file to be used into other
pipelines.

It will process the fasta input files to create indexes and dictionary.

### Panel capture

#### AlignDNAcapture

This workflow align DNA capture sequence against the reference genome.

It produce some quality metrics from different tools.

#### VariantCallingCaptureHC

This workflow apply variant calling by HaplotypeCaller on a capture sequencing.

It produce some quality metrics from different tools.

#### PanelCapture

This workflow performed an analysis for capture constitutionnal sample.
It calls the following subworkflow :

- ***PrepareGenome*** (if necessary, i.e. fasta provided without indexes)
- ***AlignDNAcapture***
- ***VariantCallingCaptureHC***

## Tools implemented

- bash :
	- findFiles
	- convertBedToIntervals
	- makeLink
	- concatenateFiles
- bgzip (v1.10)
	- bgzip (compress, decompress)
- bcftools (v1.10) :
	- index
	- merge
	- norm
	- stats
- bedtools (v2.29.2) :
	- intersect
	- sort
	- coverage
- BWA (0.7.17-r1188) :
	- mem
	- index
- Crumble (0.8.3) :
	- crumble
- Cutadapt (2.10) :
	- adaptersTrimming
- dwgsim (0.1.11) :
	- simulateReadsIllumina
- FastQC (v0.11.9) :
	- fastqc
	- fastqcNano
	- fastqcCasava
- GATK4 (4.1.8.1) :
	- ApplyBQSR
	- BaseRecalibrator
	- BedToIntervalList
	- IntervalListToBed
	- HaplotypeCaller
	- CollectMultipleMetrics
	- GatherBamFiles
	- GatherBQSRReports
	- GatherVcfFiles
	- LeftAlignIndels
	- ReorderSam
	- DepthOfCoverage (BETA)
	- SplitIntervals
	- SplitVcfs
	- VariantFiltration
	- MergeVcfs
	- SortVcf
- Sambamba (0.6.6) :
	- index
	- flagstat
	- markdup
	- sort
	- view
- samtools (1.10) :
	- bedcov
	- index
	- flagstat
	- sort
	- dict
	- view
	- faidx
	- fqidx
	- markdup
	- fixmate
	- merge
- Star :
	- genomeGenerate
- Qualimap (v.2.2.2-dev) :
	- bamqc
- tabix (1.10):
	- index
