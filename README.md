# MobiDL 2

![plateforms:Linux-64|osx-64](https://img.shields.io/badge/plateform-Linux--64|osx--64-green?&style=for-the-badge)

MobiDL 2 is a collection of tools wrapped in WDL to be used in any WDL pipelines.

## __WARNING : BETA VERSION DO NOT USE IN PRODUCTION__

## Introduction

This repo provide a set of workflow in WDL (see [workflows](#workflows-implemented))
and a set of tools wrapped in WDL tasks (see [tools](#list-of-tasks)).

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

- ***AlignDNAcapture***
- ***VariantCallingCaptureHC***

## List of tasks

<table class="tg">
<thead>
  <tr>
    <th>Tool</th>
    <th>Task</th>
    <th>Version</th>
  </tr>
</thead>
<tbody>
  <tr>
    <td rowspan="7">bash</td>
    <td>findFiles</td>
    <td>0.0.1</td>
  </tr>
  <tr>
    <td>convertBedToIntervals</td>
    <td>0.0.1</td>
  </tr>
  <tr>
    <td>makeLink</td>
    <td>0.0.1</td>
  </tr>
  <tr>
    <td>concatenateFiles</td>
    <td>0.0.1</td>
  </tr>
  <tr>
    <td>wget</td>
    <td>0.0.1</td>
  </tr>
  <tr>
    <td>gzip</td>
    <td>0.0.1</td>
  </tr>
  <tr>
    <td>sortgtf</td>
    <td>0.0.1</td>
  </tr>
  <tr>
    <td rowspan="3">bedtools</td>
    <td>intersect</td>
    <td>0.0.1</td>
  </tr>
  <tr>
    <td>sort</td>
    <td>0.0.1</td>
  </tr>
  <tr>
    <td>coverage</td>
    <td>0.0.1</td>
  </tr>
  <tr>
    <td>bgzip</td>
    <td>bgzip</td>
    <td>0.0.1</td>
  </tr>
  <tr>
    <td rowspan="4">bcftools</td>
    <td>index</td>
    <td>0.0.1</td>
  </tr>
  <tr>
    <td>merge</td>
    <td>0.0.1</td>
  </tr>
  <tr>
    <td>norm</td>
    <td>0.0.1</td>
  </tr>
  <tr>
    <td>stats</td>
    <td>0.0.1</td>
  </tr>
  <tr>
    <td rowspan="2">BWA</td>
    <td>mem</td>
    <td>0.0.1</td>
  </tr>
  <tr>
    <td>index</td>
    <td>0.0.1</td>
  </tr>
  <tr>
    <td>Crumble</td>
    <td>crumble</td>
    <td>0.0.1</td>
  </tr>
  <tr>
    <td rowspan="3">Cutadapt</td>
    <td>adaptersTrimming</td>
    <td>0.0.1</td>
  </tr>
  <tr>
    <td>qualityTrimming</td>
    <td>0.0.1</td>
  </tr>
  <tr>
    <td>hardTrimming</td>
    <td>0.0.1</td>
  </tr>
  <tr>
    <td>DWGSim</td>
    <td>simulateReadsIllumina</td>
    <td>0.0.1</td>
  </tr>
  <tr>
    <td rowspan="3">FastQC</td>
    <td>fastqc</td>
    <td>0.0.1</td>
  </tr>
  <tr>
    <td>fastqcNano</td>
    <td>0.0.1</td>
  </tr>
  <tr>
    <td>fastqcCasava</td>
    <td>0.0.1</td>
  </tr>
  <tr>
    <td rowspan="17">GATK4</td>
    <td>ApplyBQSR</td>
    <td>0.0.1</td>
  </tr>
  <tr>
    <td>BaseRecalibrator</td>
    <td>0.0.1</td>
  </tr>
  <tr>
    <td>BedToIntervalList</td>
    <td>0.0.1</td>
  </tr>
  <tr>
    <td>IntervalListToBed</td>
    <td>0.0.1</td>
  </tr>
  <tr>
    <td>HaplotypeCaller</td>
    <td>0.0.1</td>
  </tr>
  <tr>
    <td>CollectMultipleMetrics</td>
    <td>0.0.1</td>
  </tr>
  <tr>
    <td>GatherBamFiles</td>
    <td>0.0.1</td>
  </tr>
  <tr>
    <td>GatherBQSRReports</td>
    <td>0.0.1</td>
  </tr>
  <tr>
    <td>GatherVcfFiles</td>
    <td>0.0.1</td>
  </tr>
  <tr>
    <td>LeftAlignIndels</td>
    <td>0.0.1</td>
  </tr>
  <tr>
    <td>ReorderSam</td>
    <td>0.0.1</td>
  </tr>
  <tr>
    <td>DepthOfCoverage (BETA)</td>
    <td>0.0.1</td>
  </tr>
  <tr>
    <td>SplitIntervals</td>
    <td>0.0.1</td>
  </tr>
  <tr>
    <td>SplitVcfs</td>
    <td>0.0.1</td>
  </tr>
  <tr>
    <td>VariantFiltration</td>
    <td>0.0.1</td>
  </tr>
  <tr>
    <td>MergeVcfs</td>
    <td>0.0.1</td>
  </tr>
  <tr>
    <td>SortVcf</td>
    <td>0.0.1</td>
  </tr>
  <tr>
    <td>MPA</td>
    <td>mpa</td>
    <td>0.0.1</td>
  </tr>
  <tr>
    <td>Qualimap</td>
    <td>bamqc</td>
    <td>0.0.1</td>
  </tr>
  <tr>
    <td>Rsync</td>
    <td>rsync</td>
    <td>0.0.1</td>
  </tr>
  <tr>
    <td rowspan="2">Regtools</td>
    <td>junctionsExtract</td>
    <td>0.0.1</td>
  </tr>
  <tr>
    <td>junctionsAnnotate</td>
    <td>0.0.1</td>
  </tr>
  <tr>
    <td rowspan="5">Sambamba</td>
    <td>index</td>
    <td>0.0.1</td>
  </tr>
  <tr>
    <td>flagstat</td>
    <td>0.0.1</td>
  </tr>
  <tr>
    <td>markdup</td>
    <td>0.0.1</td>
  </tr>
  <tr>
    <td>sort</td>
    <td>0.0.1</td>
  </tr>
  <tr>
    <td>view</td>
    <td>0.0.1</td>
  </tr>
  <tr>
    <td rowspan="11">Samtools</td>
    <td>bedcov</td>
    <td>0.0.1</td>
  </tr>
  <tr>
    <td>index</td>
    <td>0.0.1</td>
  </tr>
  <tr>
    <td>flagstat</td>
    <td>0.0.1</td>
  </tr>
  <tr>
    <td>sort</td>
    <td>0.0.1</td>
  </tr>
  <tr>
    <td>dict</td>
    <td>0.0.1</td>
  </tr>
  <tr>
    <td>view</td>
    <td>0.0.1</td>
  </tr>
  <tr>
    <td>faidx</td>
    <td>0.0.1</td>
  </tr>
  <tr>
    <td>fqidx</td>
    <td>0.0.1</td>
  </tr>
  <tr>
    <td>markdup</td>
    <td>0.0.1</td>
  </tr>
  <tr>
    <td>fixmate</td>
    <td>0.0.1</td>
  </tr>
  <tr>
    <td>merge</td>
    <td>0.0.1</td>
  </tr>
  <tr>
    <td>Star</td>
    <td>genomeGenerate</td>
    <td>0.0.1</td>
  </tr>
  <tr>
    <td>Tabix</td>
    <td>index</td>
    <td>0.0.1</td>
  </tr>
  <tr>
    <td rowspan="3">Vardict-java</td>
    <td>vardictSoloAmplicons</td>
    <td>0.0.1</td>
  </tr>
  <tr>
    <td>teststrandbias</td>
    <td>0.0.1</td>
  </tr>
  <tr>
    <td>var2vcf_valid</td>
    <td>0.0.1</td>
  </tr>
  <tr>
    <td rowspan="2">VEP</td>
    <td>vep_cache</td>
    <td>0.0.1</td>
  </tr>
  <tr>
    <td>vep_install</td>
    <td>0.0.1</td>
  </tr>
</tbody>
</table>
