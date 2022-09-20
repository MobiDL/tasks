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
		author: "Olivier Ardouin"
		email: "o-ardouin(at)chu-montpellier.fr"
		version: "Beta.0.0.1"
		date: "2022-09-14"
	}

	input {
		String path_exe = "freebayes"

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
		set -exo pipefail
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
			description: 'Path used as executable [default: "freebayes"]',
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

task freebayes {
	meta {
		author: "Olivier Ardouin"
		email: "o-ardouin(at)chu-montpellier.fr"
		version: "Beta.0.0.1"
		date: "2022-09-19"
	}

	input {
		String path_exe = "freebayes"

		#input and output:
		File? in = ""
		Array[File]? bams
		File? bamList
		String? outputName
		String ext = ".bam"
		String? outputPath
		File refFasta
		File? Target
		Array[String]? Regions
		File? Samples
		File? Populations
		File? cnvMap
		String? Trace
		String? FailedAlleles
		File? VariantExtendedInput
		File? VariantRestrictInput
		File? HaplotypeBaseVCF
		Boolean reportAllHaplotypeAlleles = false
		Boolean reportMonomorphic = false

		#Reporting:
		Float? pvar

		#Population models:
		Float? theta
		Int? ploidy
		Boolean PooledDicrete = false
		Boolean PooledContinuous = false

		#reference allele:
		Boolean useRefAllele = false
		String? ReferenceQuality

		#allele scope:
		Boolean NoSNP = false
		Boolean NoIndels = false
		Boolean NoMNPs = false
		Boolean NoComplex = false
		Int? useBestN
		Int? HaplotypeLenght
		Int? MinRepeatSize
		Int? MinRepeatEntropy
		Boolean NoPartialObservations = false

		#indel realignment:
		Boolean DontLeftAlignIndels = false

		#input filters:
		Boolean UseDuplicateReads = false
		Int? MinMappingQuality
		Int? MinBaseQuality
		Int? MinSupportingAlleleQsum
		Int? MinSupportingMappingQsum
		Int? MismatchBaseQualityThreshold
		Int? ReadMismatchLimit
		Float? ReadMaxMismatchFraction
		Int? ReadSnpLimit
		Int? ReadIndelLimit
		Boolean standardFilters = false
		Float? MinAlternateFraction
		Int? MinAlternateCount
		Float? MinAlternateQsum
		Int? MinAlternateTotal
		Int? MinCoverage

		#population priors:
		Boolean NoPopulationPriors = false

		#mappability priors:
		Boolean hwePriorsOff = false
		Boolean binomialObsPriorsOff = false
		Boolean alleleBalancePriorsOff = false

		#genotype likelihoods:
		File? observationBias
		Float? baseQualityCap
		Float? probContamination
		Boolean legacyGLS = false
		File? contaminationEstimates

		#algorithmic features:
		Boolean reportGenotypeLikelihoodMax = false
		Int? genotypingMaxIterations
		Int? genotypingMaxBanddepth
		String? posteriorIntegrationLimits
		Boolean excludeUnobservedGenotypes = false
		Float? genotypeVariantThreshold
		Boolean useMappingQuality = false
		Boolean harmonicIndelQuality = false
		Float? readDependenceFactor
		Boolean genotypeQualities = false

		Int threads = 1
		Int memoryByThreads = 768
		String? memory
	}

	String totalMem = if defined(memory) then memory else memoryByThreads*threads + "M"
	Boolean inGiga = (sub(totalMem,"([0-9]+)(M|G)", "$2") == "G")
	Int memoryValue = sub(totalMem,"([0-9]+)(M|G)", "$1")
	Int totalMemMb = if inGiga then memoryValue*1024 else memoryValue
	Int memoryByThreadsMb = floor(totalMemMb/threads)

	# define input
	String InputBamList = if defined(bamList) then " --bam-list ~{bamList}" else ""
	String InputBamArray = if defined(bams) then "--bam " else ""

	# define output file
	String OutputName = if defined(outputName) then outputName else basename(in,ext) + ".vcf"
	String OutputFile = if defined(outputPath) then "~{outputPath}/~{outputName}" else "~{outputName}"

	#input parameters
	String InputTarget = if defined(Target) then "--targets ~{Target} " else ""
	String InputRegions = if defined(Regions) then "--region ~{Regions} " else ""
	String InputSamples = if defined(Samples) then "--samples ~{Samples} " else ""
	String InputPopulations = if defined(Populations) then "--population ~{Populations} " else ""
	String InputcnvMap = if defined(cnvMap) then "--cnv-map ~{cnvMap} " else ""
	String InputTrace = if defined(Trace) then "--trace ~{Trace} " else ""
	String InputFailedAlleles = if defined(FailedAlleles) then "--failed-alleles ~{FailedAlleles} " else ""
	String InputVariantExtendedInput = if defined(VariantExtendedInput) then "--variant-input ~{VariantExtendedInput} " else ""
	String InputVariantRestrictInput = if defined(VariantRestrictInput) then "--only-use-input-alleles ~{VariantRestrictInput} " else ""
	String InputHaplotypeBaseVCF = if defined(HaplotypeBaseVCF) then "--haplotype-basis-alleles ~{HaplotypeBaseVCF} " else ""
	String InputreportAllHaplotypeAlleles = if reportAllHaplotypeAlleles then "--report-all-haplotype-alleles ~{reportAllHaplotypeAlleles} " else ""
	String InputreportMonomorphic = if reportMonomorphic then "--report-monomorphic ~{reportMonomorphic} " else ""

	#reporting parameters
	String Inputpvar = if defined(pvar) then "--pvar ~{pvar} " else ""

	#Population model parameters
	String Inputtheta = if defined(theta) then "--theta ~{theta} " else ""
	String Inputploidy = if defined(ploidy) then "--ploidy ~{ploidy} " else ""
	String InputPooledDicrete = if PooledDicrete then "--pooled-discrete " else ""
	String InputPooledContinuous = if PooledContinuous then "--pooled-continuous " else ""

	#Reference allele parameters
	String InputuseRefAllele = if useRefAllele then "--use-reference-allele " else ""
	String InputReferenceQuality = if defined(ReferenceQuality) then "--reference-quality ~{ReferenceQuality} " else ""

	#allele scope parameters
	String InputNoSNP = if NoSNP then "--no-snps " else ""
	String InputNoIndels = if NoIndels then "--no-indels " else ""
	String InputNoMNPs = if NoMNPs then "--no-mnps " else ""
	String InputNoComplex = if NoComplex then "--no-complex " else ""
	String InputuseBestN = if defined(useBestN) then "--use-best-n-alleles ~{useBestN} " else ""
	String InputHaplotypeLenght = if defined(HaplotypeLenght) then "--haplotype-length {HaplotypeLenght} " else ""
	String InputMinRepeatSize = if defined(MinRepeatSize) then "--min-repeat-size ~{MinRepeatSize} " else ""
	String InputMinRepeatEntropy = if defined(MinRepeatEntropy) then "--min-repeat-entropy ~{MinRepeatEntropy} " else ""
	String InputNoPartialObservations = if NoPartialObservations then "--no-partial-observations " else ""

	#Indel realignement parameter
	String InputDontLeftAlignIndels = if DontLeftAlignIndels then "--dont-left-align-indels " else ""

	#input filters parameters
	String InputUseDuplicateReads = if UseDuplicateReads then "--use-duplicate-reads " else ""
	String InputMinMappingQuality = if defined(MinMappingQuality) then "--min-mapping-quality ~{MinMappingQuality} " else ""
	String InputMinBaseQuality = if defined(MinBaseQuality) then "--min-base-quality ~{MinBaseQuality} " else ""
	String InputMinSupportingAlleleQsum = if defined(MinSupportingAlleleQsum) then "--min-supporting-allele-qsum ~{MinSupportingAlleleQsum} " else ""
	String InputMinSupportingMappingQsum = if defined(MinSupportingMappingQsum) then "--min-supporting-mapping-qsum ~{MinSupportingMappingQsum} " else ""
	String InputMismatchBaseQualityThreshold = if defined(MismatchBaseQualityThreshold) then "--mismatch-base-quality-threshold ~{MismatchBaseQualityThreshold} " else ""
	String InputReadMismatchLimit = if defined(ReadMismatchLimit) then "--read-mismatch-limit ~{ReadMismatchLimit} " else ""
	String InputReadMaxMismatchFraction = if defined(ReadMaxMismatchFraction) then "--read-max-mismatch-fraction ~{ReadMaxMismatchFraction} " else ""
	String InputReadSnpLimit = if defined(ReadSnpLimit) then "--read-snp-limit ~{ReadSnpLimit} " else ""
	String InputReadIndelLimit = if defined(ReadIndelLimit) then "--read-indel-limit ~{ReadIndelLimit} " else ""
	String InputstandardFilters = if standardFilters then "--standard-filters " else ""
	String InputMinAlternateFraction = if defined(MinAlternateFraction) then "--min-alternate-fraction ~{MinAlternateFraction} " else ""
	String InputMinAlternateCount = if defined(MinAlternateCount) then "--min-alternate-count ~{MinAlternateCount} " else ""
	String InputMinAlternateQsum = if defined(MinAlternateQsum) then "--min-alternate-qsum ~{MinAlternateQsum} " else ""
	String InputMinAlternateTotal = if defined(MinAlternateTotal) then "--min-alternate-total~{MinAlternateTotal} " else ""
	String InputMinCoverage = if defined(MinCoverage) then "--min-coverage ~{MinCoverage} " else ""

	#population priors parameter
	String InputNoPopulationPriors = if NoPopulationPriors then "--no-population-priors " else ""

	#mappability priors parameters
	String InputhwePriorsOff = if hwePriorsOff then "--hwe-priors-off " else ""
	String InputbinomialObsPriorsOff = if binomialObsPriorsOff then "--binomial-obs-priors-off " else ""
	String InputalleleBalancePriorsOff = if alleleBalancePriorsOff then "--allele-balance-priors-off " else ""

	#genotype likelihoods parameters
	String InputobservationBias = if defined(observationBias) then "--observation-bias ~{observationBias} " else ""
	String InputbaseQualityCap = if defined(baseQualityCap) then "--base-quality-cap ~{baseQualityCap} " else ""
	String InputprobContamination = if defined(probContamination) then "--prob-contamination ~{probContamination} " else ""
	String InputlegacyGLS = if legacyGLS then "--legacy-gls " else ""
	String InputcontaminationEstimates = if defined(contaminationEstimates) then "--contamination-estimates ~{contaminationEstimates} " else ""

	#algorithmic features parameters
	String InputreportGenotypeLikelihoodMax = if reportGenotypeLikelihoodMax then "--report-genotype-likelihood-max " else ""
	String InputgenotypingMaxIterations = if defined(genotypingMaxIterations) then "--genotyping-max-iterations ~{genotypingMaxIterations} " else ""
	String InputgenotypingMaxBanddepth = if defined(genotypingMaxBanddepth) then "--genotyping-max-banddepth ~{genotypingMaxBanddepth} " else ""
	String InputposteriorIntegrationLimits = if defined(posteriorIntegrationLimits) then "--posterior-integration-limits ~{posteriorIntegrationLimits} " else ""
	String InputexcludeUnobservedGenotypes = if excludeUnobservedGenotypes then "--exclude-unobserved-genotypes " else ""
	String InputgenotypeVariantThreshold = if defined(genotypeVariantThreshold) then "--genotype-variant-threshold ~{genotypeVariantThreshold} " else ""
	String InputuseMappingQuality = if useMappingQuality then "--use-mapping-quality " else ""
	String InputharmonicIndelQuality = if harmonicIndelQuality then "--harmonic-indel-quality " else ""
	String InputreadDependenceFactor = if defined(readDependenceFactor) then "--read-dependence-factor ~{readDependenceFactor} " else ""
	String InputgenotypeQualities = if genotypeQualities then "--genotype-qualities " else ""

	command <<<
		set -exo pipefail
		if [[ ! -f ~{OutputFile} ]]; then
			mkdir -p $(dirname ~{OutputFile})
		fi

		if [[ -z "~{in}" && -z "~{InputBamArray}" ]]; then
			echo "Freebayes : No entry (bam) provided. exit 1";
			exit 1
		fi

		BamArray=""
		if [[ -z "~{InputBamArray}" ]]; then
			BamArray="~{InputBamArray} ~{sep=' --bam ' bams} "
		fi

		~{path_exe} " ~{in} " \
		~{BamArray} \
		~{InputBamList} \
		--vcf ~{OutputFile} \
		--fasta-reference ~{refFasta} \
		~{InputTarget} \
		~{InputRegions} \
		~{InputSamples} \
		~{InputPopulations} \
		~{InputcnvMap} \
		~{InputTrace} \
		~{InputFailedAlleles} \
		~{InputVariantExtendedInput} \
		~{VariantRestrictInput} \
		~{InputHaplotypeBaseVCF} \
		~{InputreportAllHaplotypeAlleles} \
		~{InputreportMonomorphic} \
		~{Inputpvar} \
		~{Inputtheta} \
		~{Inputploidy} \
		~{InputPooledDicrete} \
		~{InputPooledContinuous} \
		~{InputuseRefAllele} \
		~{InputReferenceQuality} \
		~{InputNoSNP} \
		~{InputNoIndels} \
		~{InputNoMNPs} \
		~{InputNoComplex} \
		~{InputuseBestN} \
		~{InputHaplotypeLenght} \
		~{InputMinRepeatSize} \
		~{InputMinRepeatEntropy} \
		~{InputNoPartialObservations} \
		~{InputDontLeftAlignIndels} \
		~{InputUseDuplicateReads} \
		~{InputMinMappingQuality} \
		~{InputMinBaseQuality} \
		~{InputMinSupportingAlleleQsum} \
		~{InputMinSupportingMappingQsum} \
		~{InputMismatchBaseQualityThreshold} \
		~{InputReadMismatchLimit} \
		~{InputReadMaxMismatchFraction} \
		~{InputReadSnpLimit} \
		~{InputReadIndelLimit} \
		~{InputstandardFilters} \
		~{InputMinAlternateFraction} \
		~{InputMinAlternateCount} \
		~{InputMinAlternateQsum} \
		~{InputMinAlternateTotal} \
		~{InputMinCoverage} \
		~{InputNoPopulationPriors} \
		~{InputhwePriorsOff} \
		~{InputbinomialObsPriorsOff} \
		~{InputalleleBalancePriorsOff} \
		~{InputobservationBias} \
		~{InputbaseQualityCap} \
		~{InputprobContamination} \
		~{InputlegacyGLS} \
		~{InputcontaminationEstimates} \
		~{InputreportGenotypeLikelihoodMax} \
		~{InputgenotypingMaxIterations} \
		~{InputgenotypingMaxBanddepth} \
		~{InputposteriorIntegrationLimits} \
		~{InputexcludeUnobservedGenotypes} \
		~{InputgenotypeVariantThreshold} \
		~{InputuseMappingQuality} \
		~{InputharmonicIndelQuality} \
		~{InputreadDependenceFactor} \
		~{InputgenotypeQualities}
	>>>

	output {
		File outputFile = "~{OutputFile}"
	}

	runtime {
		cpu: "~{threads}"
		requested_memory_mb_per_core: "${memoryByThreadsMb}"
	}

	parameter_meta {
		path_exe: {
			description: 'Path used as executable [default: "freebayes"]',
			category: 'System'
		}
		in: {
			description: 'Path of the firts (or single) bam file.',
			category: 'mandatory if no bamList'
		}
		bams: {
			description: 'Paths of the others bam file to be analyzed.',
			category: 'option'
		}
		bamList: {
			description: 'A file containing a list of BAM files to be analyzed.',
			category: 'mandatory if no in'
		}
		outputPath: {
			description: 'Output path where file (VCF) were generated.',
			category: 'Output path/name option'
		}
		refFasta: {
			description: 'Reference Sequence in FASTA format used for the analysis. An index file (refFasta.fai) will be created if none exists.',
			category: 'mandatory'
		}
		Target: {
			description: 'Limit analysis to target listed in the bed format file. If neither Target nor Region are specified, every position in RefFasta will be analyzed.',
			category: 'option'
		}
		Regions: {
			description: 'Limit analysis to the specified region. format = <chrom>:<start_position>-<end_position> , 0 base coordinate, end position not included (same as BED format). Either - or .. maybe used as a separator. If neither Target nor Region are specified, every position in RefFasta will be analyzed.',
			category: 'option'
		}
		Samples: {
			description: 'Limit analysis to samples listed (one per line) in the file. By default all samples in input BAM files will be analyzed.',
			category: 'option'
		}
		Populations: {
			description: 'Each line of FILE should list a sample and a population which it is part of. The population-based bayesian inference model will then be partitioned on the basis of the populations.',
			category: 'option'
		}
		cnvMap: {
			description: 'Read a copy number map from the BED formated file provided. The format is: reference sequence, start, end, sample name, copy number... for each region in each sample which does not have the default copy number as set by --ploidy.',
			category: 'option'
		}
		Trace: {
			descrition: 'Output an algorithmic trace to output file path provided.',
			category: 'option'
		}
		FailedAlleles: {
			description: 'Write a BED format file of the analyzed positions which do not pass --pvar to output file path provided..',
			category: 'option'
		}
		VariantExtendedInput: {
			description: 'Use variants reported in VCF formated file provided as input to the algorithm. Variants in this file will be included in the output even if there is not enough support in the data to pass input filters.',
			category: 'option'
		}
		VariantRestrictInput: {
			description: 'VCF formated file. Only provide variant calls and genotype likelihoods for sites and alleles which are provided in this VCF input, and provide output in the VCF for all input alleles, not just those which have support in the data.',
			category: 'option'
		}
		HaplotypeBaseVCF: {
			description: 'When specified, only variant alleles provided in this input VCF formated files will be used for the construction of complex or haplotype alleles.'
			category: 'option'
		}
		reportAllHaplotypeAlleles: {
			description: 'If true: At sites where genotypes are made over haplotype alleles, provide information about all alleles in output, not only those which are called.',
			category: 'option; default = false'
		}
		reportMonomorphic : {
			description: 'If true: Report even loci which appear to be monomorphic, and report all considered alleles, even those which are not in called genotypes. Loci which do not have any potential alternates have . for ALT.',
			category: 'option'
		}
		pvar: {
			description: 'Report sites if the probability that there is a polymorphism at the site is greater than value provided. Default is 0.0. Note that post- filtering is generally recommended over the use of this parameter.',
			category: 'option'
		}
		theta: {
			description: 'The expected mutation rate or pairwise nucleotide diversity among the population under analysis. This serves as the single parameter to the Ewens Sampling Formula prior model default: 0.001',
			category: 'option'
		}
		ploidy: {
			description: 'Sets the default ploidy for the analysis to N. Default = 2',
			category: 'option'
		}
		PooledDicrete: {
			description: 'If true : Assume that samples result from pooled sequencing. Model pooled samples using discrete genotypes across pools. When using this flag, set --ploidy to the number of alleles in each sample or use the --cnv-map to define per-sample ploidy. Default is false',
			category: 'option'
		}
		PooledContinuous: {
			description: 'If true : Output all alleles which pass input filters, regardles of genotyping outcome or model. Default = false.',
			category: 'option'
		}
		useRefAllele: {
			description: 'If true includes the reference allele in the analysis as if it is another sample from the same population. Default = false.',
			category: 'option'
		}
		ReferenceQuality: {
			description: 'MQ,BQ : Assign mapping quality of MQ to the reference allele at each : site and base quality of BQ. Default is 100,60',
			category: 'option'
		}
		NoSNP: {
			description: 'If true : Ignore SNP alleles. Default is false',
			category: 'option'
		}
		NoIndels: {
			description: 'If true : Ignore insertion and deletion alleles. Default is false',
			category: 'option'
		}
		NoMNPs: {
			description: 'If true : Ignore multi-nuceotide polymorphisms, MNPs. Default is false',
			category: 'option'
		}
		NoComplex: {
			description: 'If true : Ignore complex events (composites of other classes). Default is false',
			category: 'option'
		}
		useBestN: {
			description: 'Evaluate only the best N SNP alleles, ranked by sum of supporting quality scores. (Set to 0 to use all; default is all).',
			category: 'option'
		}
		HaplotypeLenght: {
			description: 'Allow haplotype calls with contiguous embedded matches of up to this length. default is 3.',
			category: 'option'
		}
		MinRepeatSize: {
			description: 'When assembling observations across repeats, require the total repeat length at least this many bp. Default is 5.',
			category: 'option'
		}
		MinRepeatEntropy:{
			description: 'To detect interrupted repeats, build across sequence until it has entropy > N bits per bp. Default is 0 : off.',
			category: 'option'
		}
		NoPartialObservations: {
			description: 'If true : Exclude observations which do not fully span the dynamically-determined detection window. Default is false: use all observations, dividing partial support across matching haplotypes when generating haplotypes.',
			category: 'option'
		}
		DontLeftAlignIndels: {
			description: 'If true : Turn off left-alignment of indels, which is enabled (false) by default.',
			category: 'option'
		}
		UseDuplicateReads: {
			description: 'If true : Include duplicate-marked alignments in the analysis. Default is false : exclude duplicates marked as such in alignments.',
			category: 'option'
		}
		MinMappingQuality:{
			description: 'Exclude alignments from analysis if they have a mapping quality less than provided value. Default is 1.',
		}
		MinBaseQuality: {
			description: 'Exclude alleles from analysis if their supporting base quality is less than provided value. Default is 0.',
			category: 'option'
		}
		MinSupportingAlleleQsum: {
			description: 'Consider any allele in which the sum of qualities of supporting observations is at least provided value. Default is 0.',
			category: 'option'
		}
		MinSupportingMappingQsum: {
			description: 'Consider any allele in which and the sum of mapping qualities of supporting reads is at least provided value. Default is 0.',
			category: 'option'
		}
		MismatchBaseQualityThreshold: {
			description: 'Count mismatches toward ReadMismatchLimit option if the base quality of the mismatch is >= provided value. Default is 10.',
			category: 'option'
		}
		ReadMismatchLimit: {
			description: 'Exclude reads with more than N (provided value) mismatches where each mismatch has base quality >= MismatchBaseQualityThreshold. Default is ~unbounded.',
			category: 'option'
		}
		ReadMaxMismatchFraction: {
			description: 'Exclude reads with more than provided value [0,1] fraction of mismatches where each mismatch has base quality >= MismatchBaseQualityThreshold. Default is 1.0 .',
			category: 'option'
		}
		ReadSnpLimit: {
			description: 'Exclude reads with more than N (provided value) base mismatches, ignoring gaps with quality >= MismatchBaseQualityThreshold Default is ~unbounded',
			category: 'option'
		}
		ReadIndelLimit: {
			description: 'Exclude reads with more than N (provided value) separate gaps. Default is ~unbounded.',
			category: 'option'
		}
		standardFilters: {
			description: 'If true : Use stringent input base and mapping quality filters Equivalent to MinMappingQuality = 30; MinBaseQuality = 20; MinSupportingAlleleQsum = 0; genotypeVariantThreshold = 0',
			category: 'option'
		}
		MinAlternateFraction: {
			description: 'Require at least this fraction of observations supporting an alternate allele within a single individual in the in order to evaluate the position. Default is 0.2 .',
			category: 'option'
		}
		MinAlternateCount: {
			description: 'Require at least this count of observations supporting an alternate allele within a single individual in order to evaluate the position. Default is 2 .',
			category: 'option'
		}
		MinAlternateQsum: {
			description: 'Require at least this sum of quality of observations supporting an alternate allele within a single individual in order to evaluate the position. Default is 0 .',
			category: 'option'
		}
		MinAlternateTotal: {
			description: 'Require at least this count of observations supporting an alternate allele within the total population in order to use the allele in analysis. Default is 1 .',
			category: 'option'
		}
		MinCoverage: {
			description: 'Require at least this coverage to process a site. Default is 0 .',
			category: 'option'
		}
		NoPopulationPriors: {
			description: 'If true : Equivalent to PooledDicrete = true ; hwePriorsOff and removal of Ewens Sampling Formula component of priors. Default is false.',
			category: 'option'
		}
		hwePriorsOff: {
			description: 'If true : Disable estimation of the probability of the combination arising under HWE given the allele frequency as estimated by observation frequency. Default is false.',
			category: 'option'
		}
		binomialObsPriorsOff: {
			description: 'If true : Disable incorporation of prior expectations about observations. Uses read placement probability, strand balance probability, and read position (5\'-3\') probability. Default is false.',
			category: 'option'
		}
		alleleBalancePriorsOff: {
			description: 'If true : Disable use of aggregate probability of observation balance between alleles as a component of the priors. Default is false.',
			category: 'option'
		}
		observationBias: {
			description: 'Read length-dependent allele observation biases from provided file. The format is [length] [alignment efficiency relative to referprovided file] where the efficiency is 1 if there is no relative observation bias.',
			category: 'option'
		}
		baseQualityCap: {
			description: 'Limit estimated observation quality by capping base quality at provided value.',
			category: 'option'
		}
		probContamination: {
			description: 'An estimate of contamination to use for all samples. Default value is 10e-9',
			category: 'option'
		}
		legacyGLS: {
			description: 'If true : Use legacy (polybayes equivalent) genotype likelihood calculations. Default is false.',
			category: 'option'
		}
		contaminationEstimates: {
			description: 'A file containing per-sample estimates of contamination, such as those generated by VerifyBamID. The format should be: sample p(read=R|genotype=AR) p(read=A|genotype=AA) Sample \'*\' can be used to set default contamination estimates.',
			category: 'option'
		}
		reportGenotypeLikelihoodMax: {
			description: 'If true : Report genotypes using the maximum-likelihood estimate provided from genotype likelihoods. Default is false.',
			category: 'option'
		}
		genotypingMaxIterations: {
			description: 'Iterate no more than N times during genotyping step. Default is 1000.',
			category: 'option'
		}
		genotypingMaxBanddepth: {
			description: 'Integrate no deeper than the Nth best genotype by likelihood when genotyping. Default is 6.',
			category: 'option'
		}
		posteriorIntegrationLimits: {
			description: 'N,M : Integrate all genotype combinations in our posterior space which include no more than N samples with their Mth best data likelihood. Default is 1,3.',
			category: 'option'
		}
		excludeUnobservedGenotypes: {
			description: 'If true : Skip sample genotypings for which the sample has no supporting reads. Default is false.',
			category: 'option'
		}
		genotypeVariantThreshold: {
			description: 'Limit posterior integration to samples where the second-best genotype likelihood is no more than log(provided value) from the highest genotype likelihood for the sample. Default is ~unbounded.',
			category: 'option'
		}
		useMappingQuality: {
			description: 'If true : Use mapping quality of alleles when calculating data likelihoods. Default is false.',
			category: 'option'
		}
		harmonicIndelQuality: {
			description: 'If true : Use a weighted sum of base qualities around an indel, scaled by the distance from the indel. By default use a minimum BQ in flanking sequence. Default is false.',
			category: 'option'
		}
		readDependenceFactor: {
			description: 'Incorporate non-independence of reads by scaling successive observations by this factor during data likelihood calculations. Default is 0.9',
			category: 'option'
		}
		genotypeQualities: {
			description: 'Calculate the marginal probability of genotypes and report as GQ in each sample field in the VCF output.',
			category: 'option'
		}
	}
}
