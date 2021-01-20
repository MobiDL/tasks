# Changelog

All notable changes to this project will be documented in this file.

## [0.0.1] - dev

### Added

#### Workflows

- Preprocess :
  - getAnnotationsFiles
  - getGenomeDir
  - preprocess
- Subworkflows :
  - alignmentShortReadsDNA
  - postProcessVCF

#### Tasks

- bash :
	- findFiles
	- convertBedToIntervals
	- makeLink
	- concatenateFiles
	- wget
	- gzip (compress/decompress)
	- sortgtf
	- fai2bed
- bcftools :
	- index
	- merge
	- norm
	- stats
- bedtools :
	- intersect
	- sort
	- coverage
- bgzip :
	- compress
	- decompress
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
- MPA :
	- mpa
- Qualimap :
	- bamqc
- regtools :
	- junctionsExtract
	- junctionsAnnotate
- rsync :
	- rsync
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
- snpEff :
	- install
- Star :
	- genomeGenerate
- tabix :
	- index
- vardict-java :
	- vardictSoloAmplicons
	- teststrandbias
	- var2vcf_valid
- vep :
	- vep_cache
	- install
- regtools :
	- junctionsExtract
	- junctionsAnnotate
