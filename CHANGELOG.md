# Changelog

All notable changes to this project will be documented in this file.

## [0.0.1] - dev

### Added

#### Workflows

- prepareGenome

#### Tasks

- bash :
	- findFiles
	- convertBedToIntervals
	- makeLink
- bgzip :
	- bgzip (compress and decompress)
- bcftools :
	- index
	- merge
	- norm
- bedtools :
	- intersect
	- sort
	- coverage
- bwa :
	- mem
	- index
- crumble :
	- crumble
- dwgsim :
	- simulateReadsIllumina
- FastQC :
	- fastqc
	- fastqcNano
	- fastqcCasava
- GATK4 :
	- ApplyBQSR
	- BaseRecalibrator
	- BedToIntervalList
	- IntervalListToBed
	- HaplotypeCaller
	- CollectMultipleMetrics
	- GatherBamFiles
	- GatherBQSRReports
	- LeftAlignIndels
	- ReorderSam
	- DepthOfCoverage (BETA)
	- SplitIntervals
- Sambamba :
	- index
	- flagstat
	- markdup
	- sort
- Samtools :
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
- Star :
	- genomeGenerate
- Qualimap :
	- bamqc
- tabix :
	- index
