# Workflow : ONT - VariantCalling

This workflow do VariantCalling on Oxford Nanopore Technologies data.

## Installation

You will need to create a new environment based on conda.

```bash
conda env create -f environment.yml
conda activate ONT-VariantCalling
```

Following specifications on [Clair](https://github.com/HKU-BAL/Clair#option-1-bioconda),
you will need to install `intervaltree` package :

```bash
pypy3 -m ensurepip
pypy3 -m pip install --no-cache-dir intervaltree==3.0.2
conda deactivate
```

## Configure your inputs

You can create your input file replacing editing the template or creating your own inputs file.

### Minimal

```json
{
	"variantCallingONT.fastqPath": "String",
	"variantCallingONT.refFa": "File",
	"variantCallingONT.refFai": "File",
	"variantCallingONT.modelPath": "String",
	"variantCallingONT.outputPath": "String",
}
```

### Extended

A full option templates is provided (inputs.json.tpl).

This template is separating in 4 categories (blank line) :
1. Global pipeline inputs (i.e. minimal)
2. Global pipeline options
3. Specific tasks inputs
4. Specific tasks options

## Launch

### Local

```bash
PATH_MOBIDL2="/path/to/MobiDL2/"
conda activate ONT-VariantCalling
cromwell run \
	-Dconfig.file=${PATH_MOBIDL2}/backends.conf/local.conf \
	--inputs /path/to/inputs.json \
	${PATH_MOBIDL2}/workflows/ONT/VariantCalling.wdl
```

### Cluster

- [TODO]

### Backends configuration

	MobiDL2 provide minimal backends configuration file.
	You can edit them to adapt to match your computing ressources.
