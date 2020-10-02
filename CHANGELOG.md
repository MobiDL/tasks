# Changelog

All notable changes to this project will be documented in this file.

## [0.0.1] - dev

### Added

#### Workflows

- prepareGenome
- alignDNAcapture
- variantCallingCaptureHC
- panelCapture

#### Tasks

- bash :
	- findFiles
	- convertBedToIntervals
	- makeLink
	- concatenateFiles
- bgzip :
	- bgzip (compress and decompress)
- bcftools :
	- index
	- merge
	- norm
	- stats
- bedtools :
	- intersect
	- sort
	- coverage
- bwa :
	- mem
	- index
- crumble :
	- crumble
- cutadapt :
	- adaptersTrimming
	- qualityTrimming
	- hardTrimming
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
	- GatherVcfFiles
	- LeftAlignIndels
	- ReorderSam
	- DepthOfCoverage (BETA)
	- SplitIntervals
	- SplitVcfs
	- VariantFiltration
	- MergeVcfs
	- SortVcf
- Sambamba :
	- index
	- flagstat
	- markdup
	- sort
	- view
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
	- merge
- Star :
	- genomeGenerate
- Qualimap :
	- bamqc
- tabix :
	- index
- vardict-java :
	- vardictSoloAmplicons
	- teststrandbias
	- var2vcf_valid
- vep :
	- vep_cache
