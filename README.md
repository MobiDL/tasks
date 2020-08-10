# MobiDL 2

MobiDL 2 is a collection of tools wrapped in WDL to be used in any WDL pipelines.

## __WARNING : BETA VERSION DO NOT USE IN PRODUCTION__

## Introduction

This repo provide a set of tools wrapped in WDL tasks (see [tools](#tools-implemented)).

The code follows WDL specifications as much as possible ([WDL-spec](https://github.com/openwdl/wdl/blob/main/versions/1.0/SPEC.md)).

## Todo list

- [ ] Make a test workflow
- [ ] Define structures for referenceGenome
	- [ ] Change in all task to use this structures
- [ ] Implement options for cluster

## Tools implemented

- bash :
	- findFiles
	- bgzip (compress, decompress)
	- convertBedToIntervals
- bcftools (v1.9) :
	- index
	- merge
	- norm
- bedtools (v2.29.2) :
	- intersect
	- sort
	- coverage
- BWA (0.7.17-r1188) :
	- mem
	- index
- Crumble (0.8.3) :
	- crumble
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
	- CollectMultipleMetrics
	- GatherBamFiles
	- GatherBQSRReports
	- LeftAlignIndels
	- ReorderSam
	- DepthOfCoverage (BETA)
	- SplitIntervals
- Sambamba (0.6.6) :
	- index
	- flagstat
	- markdup
	- sort
- samtools (1.9) :
	- bedcov
	- index
	- flagstat
	- sort
	- dict
	- view
	- faidx
	- fqidx
- Qualimap (v.2.2.2-dev) :
	- bamqc
