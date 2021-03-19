version 1.0

# MobiDL 2.0 - MobiDL 2 is a collection of tools wrapped in WDL to be used in any WDL pipelines.
# Copyright (C) 2021 MoBiDiC
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <https://www.gnu.org/licenses/>.

task get_version {
	meta {
		author: "David BAUX"
		email: "d-baux(at)chu-montpellier.fr"
		version: "0.0.1"
		date: "2021-03-19"
	}

	input {
		String path_exe = "longshot"

		Int threads = 1
		Int memoryByThreads = 768
		String? memory
	}

	String totalMem = if defined(memory) then memory else memoryByThreads*threads + "M"
	Boolean inGiga = (sub(totalMem,"([0-9]+)(M|G)", "$2") == "G")
	Int memoryValue = sub(totalMem,"([0-9]+)(M|G)", "$1")
	Int totalMemMb = if inGiga then memoryValue*1024 else memoryValue
	Int memoryByThreadsMb = floor(totalMemMb/threads)

	command <<<
		~{path_exe} --version
	>>>

	output {
		String version = read_string(stdout())
	}

	runtime {
		cpu: "~{threads}"
		requested_memory_mb_per_core: "${memoryByThreadsMb}"
	}

	parameter_meta {
		path_exe: {
			description: 'Path used as executable [default: "mytool"]',
			category: 'System'
		}
		threads: {
			description: 'Sets the number of threads [default: 1]',
			category: 'System'
		}
		memory: {
			description: 'Sets the total memory to use ; with suffix M/G [default: (memoryByThreads*threads)M]',
			category: 'System'
		}
		memoryByThreads: {
			description: 'Sets the total memory to use (in M) [default: 768]',
			category: 'System'
		}
	}
}

task module1 {
	meta {
		author: "Charles VAN GOETHEM"
		email: "c-vangoethem(at)chu-montpellier.fr"
		version: "0.0.1"
		date: "2021-03-19"
	}

	input {
		String path_exe = "mytool"

		File input

		String? outputPath
		String? name
		String subString = "(.+\.bam)"
		String subStringReplace = ""

    File bamFile
    File bamFileIndex
    File refGenome
    File refGenomeIndex
    String? region
    String outputVCF
    File potentialVariants
    String outBam
    Int minCov = 6
    Int maxCov = 8000
    Int minMapq = 20
    Float minAlleleQual = 7.0
    Float hapAssignmentQual = 20.0
    Float potentialSNVCutoff = 20.0
    Int minAltCount = 3
    Float minAltFrac = 0.125
    Float hapConvergeDelta = 0.0001
    Int anchorLength = 6
    Int maxSNVs = 3
    Int maxWindow = 50
    Int maxCigarIndel = 20
    Int bandWidth = 20
    String densityParams = "10:500:50"
    String? sampleId
    Float homSnvRate = 0.0005
    Float hetSnvRate = 0.001
    Float tsTvRatio = 0.5
    Float strandBiasPvalueCutoff = 0.01
    String? variantDebugDir

		Int threads = 1
		Int memoryByThreads = 768
		String? memory
	}

	String totalMem = if defined(memory) then memory else memoryByThreads*threads + "M"
	Boolean inGiga = (sub(totalMem,"([0-9]+)(M|G)", "$2") == "G")
	Int memoryValue = sub(totalMem,"([0-9]+)(M|G)", "$1")
	Int totalMemMb = if inGiga then memoryValue*1024 else memoryValue
	Int memoryByThreadsMb = floor(totalMemMb/threads)

	String baseName = if defined(name) then name else sub(basename(fastqR1),subString,subStringReplace)
	String outputFile = if defined(outputPath) then "~{outputPath}/~{baseName}.bam" else "~{baseName}.bam"

	command <<<
		~{path_exe} \
    --bam ~{bamFile} \
    --ref ~{refGenome} \
    ~{default="" "--region " + region} \
    --out ~{outputVCF} \
    ~{default="" "--potential_variants" + potentialVariants} \
    ~{default="" "--out_bam " + outBam} \
    ~{default="" "--min_cov " + minCov} \
    ~{default="" "--max_cov " + maxCov} \
    ~{default="" "--min_mapq" + minMapq} \
    ~{default="" "--min_allele_qual " + minAlleleQual} \
    ~{default="" "--hap_assignment_qual " + hapAssignmentQual} \
    ~{default="" "--potential_snv_cutoff " + potentialSNVCutoff} \
    ~{default="" "--min_alt_count " + minAltCount} \
    ~{default="" "--min_alt_frac " + minAltFrac} \
    ~{default="" "--hap_converge_delta " + hapConvergeDelta} \
    ~{default="" "--anchor_length " + anchorLength} \
    ~{default="" "--max_snvs " + maxSNVs} \
    ~{default="" "--max_window " + maxWindow} \
    ~{default="" "--max_cigar_indel " + maxCigarIndel} \
    ~{default="" "--band_width " + bandWidth} \
    ~{default="" "--density_params " + densityParams} \
    ~{default="" "--sample_id " + sampleId} \
    ~{default="" "--hom_snv_rate " + homSnvRate} \
    ~{default="" "--het_snv_rate " + hetSnvRate} \
    ~{default="" "--ts_tv_ratio " + tsTvRatio} \
    ~{default="" "--strand_bias_pvalue_cutoff " + strandBiasPvalueCutoff} \
    ~{default="" "--variant_debug_dir " + variantDebugDir} \
	>>>

	output {
		File output = "~{outputFile}"
	}

	runtime {
		cpu: "~{threads}"
		requested_memory_mb_per_core: "${memoryByThreadsMb}"
	}

	# generate documentation
	parameter_meta {
		path_exe: {
			description: 'Path used as executable [default: "mytool"]',
			category: 'System'
		}
		threads: {
			description: 'Sets the number of threads [default: 1]',
			category: 'System'
		}
		memory: {
			description: 'Sets the total memory to use ; with suffix M/G [default: (memoryByThreads*threads)M]',
			category: 'System'
		}
		memoryByThreads: {
			description: 'Sets the total memory to use (in M) [default: 768]',
			category: 'System'
		}
    bamFile: {
      description: 'sorted, indexed BAM file with error-prone reads',
      category: 'input'
    }
    refGenome: {
      description: 'indexed FASTA reference that BAM file is aligned to',
      category: 'input'
    }
    region: {
      description: 'Region in format <chrom> or <chrom:start-stop> in which to call variants (1-based, inclusive)',
      category: 'input'
    }
    outputVCF: {
      description: 'output VCF file with called variants',
      category: 'output'
    }
    potentialVariants: {
      description: 'Genotype and phase the variants in this tabixed/bgzipped VCF instead of using pileup method to find variants. Triallelic variants and structural variants are currently not supported.',
      category: 'input'
    }
    outBam: {
      description: 'Write new bam file with haplotype tags (HP:i:1 and HP:i:2) for reads assigned to each haplotype, any existing HP and PS tags are removed',
      category: 'option'
    }
    minCov: {
      description: 'Minimum coverage (of reads passing filters) to consider position as a potential SNV. [default: 6]',
      category: 'option'
    }
    maxCov: {
      description: 'Maximum coverage (of reads passing filters) to consider position as a potential SNV. [default: 8000]',
      category: 'option'
    }
    minMapq: {
      description: 'Minimum mapping quality to use a read. [default: 20]',
      category: 'option'
    }
    minAlleleQual: {
      description: 'Minimum estimated quality (Phred-scaled) of allele observation on read to use for genotyping/haplotyping. [default: 7.0]',
      category: 'option'
    }
    hapAssignmentQual: {
      description: 'Minimum quality (Phred-scaled) of read->haplotype assignment (for read separation). [default: 20.0]',
      category: 'option'
    }
    potentialSNVCutoff: {
      description: 'Consider a site as a potential SNV if the original PHRED-scaled QUAL score for 0/0 genotype is below this amount (a larger value considers more potential SNV sites). [default: 20.0]',
      category: 'option'
    }
    minAltCount: {
      description: 'Require a potential SNV to have at least this many alternate allele observations. [default: 3]',
      category: 'option'
    }
    minAltFrac: {
      description: 'Require a potential SNV to have at least this fraction of alternate allele observations. [default: 0.125]',
      category: 'option'
    }
    hapConvergeDelta: {
      description: 'Terminate the haplotype/genotype iteration when the relative change in log-likelihood falls below this amount. Setting a larger value results in faster termination but potentially less accurate results. [default:0.0001]',
      category: 'option'
    }
    anchorLength: {
      description: 'Length of indel-free anchor sequence on the left and right side of read realignment window. [default: 6]',
      category: 'option'
    }
    maxSNVs: {
      description: 'Cut off variant clusters after this many variants. 2^m haplotypes must be aligned against per read for a variant cluster of size m. [default: 3]',
      category: 'option'
    }
    maxWindow: {
      description: 'Maximum "padding" bases on either side of variant realignment window [default: 50]',
      category: 'option'
    }
    maxCigarIndel: {
      description: 'Throw away a read-variant during allelotyping if there is a CIGAR indel (I/D/N) longer than this amount in its window. [default: 20]',
      category: 'option'
    }
    bandWidth: {
      description: 'Minimum width of alignment band. Band will increase in size if sequences are different lengths. [default: 20]',
      category: 'option'
    }
    densityParams: {
      description: 'Parameters to flag a variant as part of a "dense cluster". Format <n>:<l>:<gq>. If there are at least n variants within l base pairs with genotype quality >=gq, then these variants are flagged as "dn" [default:10:500:50]',
      category: 'option'
    }
    sampleId: {
      description: 'Specify a sample ID to write to the output VCF [default: SAMPLE]',
      category: 'option'
    }
    homSnvRate: {
      description: 'Specify the homozygous SNV Rate for genotype prior estimation [default:0.0005]',
      category: 'option'
    }
    hetSnvRate: {
      description: 'Specify the heterozygous SNV Rate for genotype prior estimation [default:0.001]',
      category: 'option'
    }
    tsTvRatio: {
      description: 'Specify the transition/transversion rate for genotype grior estimation [default: 0.5]',
      category: 'option'
    }
    strandBiasPvalueCutoff: {
      description: 'Remove a variant if the allele observations are biased toward one strand (forward or reverse) according to Fisher\'s exact test. Use this cutoff for the two-tailed P-value. [default: 0.01]',
      category: 'option'
    }
    variantDebugDir: {
      description: 'write out current information about variants at each step of algorithm to files in this directory',
      category: 'option'
    }
	}
}
