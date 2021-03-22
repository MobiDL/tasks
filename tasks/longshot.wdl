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

task longshot {
	meta {
		author: "David BAUX"
		email: "d-baux(at)chu-montpellier.fr"
		version: "0.0.2"
		date: "2021-03-22"
	}

	input {
		String path_exe = "longshot"

		File input

		String? outputPath
		String? name
		String subString = "(\.bam)"
		String subStringReplace = ""

		Boolean autoMaxCov = false
		Boolean stableAlignment = false
		Boolean forceOverwrite = false
		Boolean maxAlignment = false
		Boolean noHaps = false
		Boolean outputRef = false

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

		if [[ ! -d $(dirname ~{outputFile}) ]]; then
			mkdir -p $(dirname ~{outputFile})
		fi

		~{path_exe} \
			--bam ~{bamFile} \
			--ref ~{refGenome} \
			~{default="" "--region " + region} \
			--out ~{outputVCF} \
			~{true="--auto_max_cov" false="" autoMaxCov} \
			~{true="--stable_alignment" false="" stableAlignment} \
			~{true="--force_overwrite" false="" forceOverwrite} \
			~{true="--max_alignment" false="" maxAlignment} \
			~{true="--no_haps" false="" noHaps} \
			~{true="--output-ref" false="" outputRef} \
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
		autoMaxCov: {
			ecription: 'Automatically calculate mean coverage for region and set max coverage to mean_coverage + 5*sqrt(mean_coverage). (SLOWER)',
			ategory : 'flag'
		}
		stableAlignment: {
			ecription: 'Use numerically-stable (logspace) pair HMM forward algorithm. Is significantly slower but may be more accurate. Tests have shown this not to be necessary for highly error prone reads (PacBio CLR).',
			ategory : 'flag'
		}
		forceOverwrite: {
			ecription: 'If output files (VCF or variant debug directory) exist, delete and overwrite them.',
			ategory : 'flag'
		}
		maxAlignment: {
			ecription: 'Use max scoring alignment algorithm rather than pair HMM forward algorithm.',
			ategory : 'flag'
		}
		noHaps: {
			ecription: 'Don\'t call HapCUT2 to phase variants.',
			ategory : 'flag'
		}
		outputRef: {
			ecription: 'print reference genotypes (non-variant), use this option only in combination with -v option.',
			ategory : 'flag'
		}
		bamFile: {
			escription: 'sorted, indexed BAM file with error-prone reads',
			ategory: 'input'
		}
		refGenome: {
			escription: 'indexed FASTA reference that BAM file is aligned to',
			ategory: 'input'
		}
		region: {
			escription: 'Region in format <chrom> or <chrom:start-stop> in which to call variants (1-based, inclusive)',
			ategory: 'input'
		}
		outputVCF: {
			escription: 'output VCF file with called variants',
			ategory: 'output'
		}
		potentialVariants: {
			escription: 'Genotype and phase the variants in this tabixed/bgzipped VCF instead of using pileup method to find variants. Triallelic variants and structural variants are currently not supported.',
			ategory: 'input'
		}
		outBam: {
			escription: 'Write new bam file with haplotype tags (HP:i:1 and HP:i:2) for reads assigned to each haplotype, any existing HP and PS tags are removed',
			ategory: 'option'
		}
		minCov: {
			escription: 'Minimum coverage (of reads passing filters) to consider position as a potential SNV. [default: 6]',
			ategory: 'option'
		}
		maxCov: {
			escription: 'Maximum coverage (of reads passing filters) to consider position as a potential SNV. [default: 8000]',
			ategory: 'option'
		}
		minMapq: {
			escription: 'Minimum mapping quality to use a read. [default: 20]',
			ategory: 'option'
		}
		minAlleleQual: {
			escription: 'Minimum estimated quality (Phred-scaled) of allele observation on read to use for genotyping/haplotyping. [default: 7.0]',
			ategory: 'option'
		}
		hapAssignmentQual: {
			escription: 'Minimum quality (Phred-scaled) of read->haplotype assignment (for read separation). [default: 20.0]',
			ategory: 'option'
		}
		potentialSNVCutoff: {
			escription: 'Consider a site as a potential SNV if the original PHRED-scaled QUAL score for 0/0 genotype is below this amount (a larger value considers more potential SNV sites). [default: 20.0]',
			ategory: 'option'
		}
		minAltCount: {
			escription: 'Require a potential SNV to have at least this many alternate allele observations. [default: 3]',
			ategory: 'option'
		}
		minAltFrac: {
			escription: 'Require a potential SNV to have at least this fraction of alternate allele observations. [default: 0.125]',
			ategory: 'option'
		}
		hapConvergeDelta: {
			escription: 'Terminate the haplotype/genotype iteration when the relative change in log-likelihood falls below this amount. Setting a larger value results in faster termination but potentially less accurate results. [default:0.0001]',
			ategory: 'option'
		}
		anchorLength: {
			escription: 'Length of indel-free anchor sequence on the left and right side of read realignment window. [default: 6]',
			ategory: 'option'
		}
		maxSNVs: {
			escription: 'Cut off variant clusters after this many variants. 2^m haplotypes must be aligned against per read for a variant cluster of size m. [default: 3]',
			ategory: 'option'
		}
		maxWindow: {
			escription: 'Maximum "padding" bases on either side of variant realignment window [default: 50]',
			ategory: 'option'
		}
		maxCigarIndel: {
			escription: 'Throw away a read-variant during allelotyping if there is a CIGAR indel (I/D/N) longer than this amount in its window. [default: 20]',
			ategory: 'option'
		}
		bandWidth: {
			escription: 'Minimum width of alignment band. Band will increase in size if sequences are different lengths. [default: 20]',
			ategory: 'option'
		}
		densityParams: {
			escription: 'Parameters to flag a variant as part of a "dense cluster". Format <n>:<l>:<gq>. If there are at least n variants within l base pairs with genotype quality >=gq, then these variants are flagged as "dn" [default:10:500:50]',
			ategory: 'option'
		}
		sampleId: {
			escription: 'Specify a sample ID to write to the output VCF [default: SAMPLE]',
			ategory: 'option'
		}
		homSnvRate: {
			escription: 'Specify the homozygous SNV Rate for genotype prior estimation [default:0.0005]',
			ategory: 'option'
		}
		hetSnvRate: {
			escription: 'Specify the heterozygous SNV Rate for genotype prior estimation [default:0.001]',
			ategory: 'option'
		}
		tsTvRatio: {
			escription: 'Specify the transition/transversion rate for genotype grior estimation [default: 0.5]',
			ategory: 'option'
		}
		strandBiasPvalueCutoff: {
			escription: 'Remove a variant if the allele observations are biased toward one strand (forward or reverse) according to Fisher\'s exact test. Use this cutoff for the two-tailed P-value. [default: 0.01]',
			ategory: 'option'
		}
		variantDebugDir: {
			escription: 'write out current information about variants at each step of algorithm to files in this directory',
			ategory: 'option'
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
