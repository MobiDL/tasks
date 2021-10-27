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


task achab {
  meta {
		author: "Olivier Ardouin"
		email: "o-ardouin(at)chu-montpellier.fr"
		version: "0.0.1"
		date: "2021-10-26"
	}
  input {
    ## sample spécific
    String SampleID
    String OutDir
    String? CaseSample
    String? FatherSample
    String? MotherSample
    String? Affected
    File mpavcf
    File? PhenolyzerFile
    String? GenesOfInterest
    Float AllelicFrequency = 0.01
    String FilterList = "PASS"
    String? CnvGeneList
    String? CustomVCF
    Float MozaicRate = 0.2
    Int MozaicDP = 5
    Boolean NewHope = false
    String? CustomInfo
    String? FavouriteGeneRef
    String? FilterCustomVCF
    String? FilterCustomVCFRegex
    Array[String]? PooledSamples
    String? SampleSubset
    Boolean AddCaseDepth = false
    File? IntersectVCF
    File? PoorCoverageFile
    File? Genemap2File
    Boolean skipCaseWT = false

    ## sytem spécific
    File AchabExe = "/mnt/Bioinfo/Softs/src/Captain-ACHAB/achab.pl"
    String PerlExe = "perl"
###    String? IDSNP --IDSNP "${IDSNP}" \  <--- missing ????
    ## run time
    Int threads = 1
		Int memoryByThreads = 768
		String? memory
  }

  String Case = if defined(CaseSample) then CaseSample else SampleID
  String Dad = if defined(FatherSample) then "--dad \"~{FatherSample}\" " else ""
  String Mum = if defined(MotherSample) then "--mum \"~{MotherSample}\" " else ""
  String affected = if defined(Affected) then "--affected ~{Affected} " else ""
  String candidates = if defined(GenesOfInterest) then "--candidates ~{GenesOfInterest} " else ""
  String cngGL = if defined(CnvGeneList) then "--cnvGeneList ~{CnvGeneList} " else ""
  String hope = if NewHope then "--newHope " else ""
  String customVcf = if defined(CustomVCF) then "--customVCF ~{CustomVCF} " else ""
  String customInfoList = if defined(CustomInfo) then "--customInfoList ~{CustomInfo} " else ""
  String favGenRef = if defined(FavouriteGeneRef) then "--favouriteGeneRef ~{FavouriteGeneRef} " else ""
  String filtCustVcf = if defined(FilterCustomVCF) then "--filterCustomVCF ~{FilterCustomVCF} " else ""
  String filtCustVcfReg = if defined(FilterCustomVCFRegex) then "--filterCustomVCFRegex ~{FilterCustomVCFRegex} " else ""
  String poolSample = if defined(PooledSamples) then "--pooledSamples " else ""
  String sampSub = if defined(SampleSubset) then "--sampleSubset ~{SampleSubset} " else ""
  String addCasDep = if AddCaseDepth then "--addCaseDepth " else ""
  String interVcf = if defined(IntersectVCF) then "--intersectVCF ~{IntersectVCF} " else ""
  String poorCov = if (defined(PoorCoverageFile) && defined(Genemap2File)) then "--poorCoverageFile ~{PoorCoverageFile} --genemap2File ~{Genemap2File} " else ""
  String SkipCase = if skipCaseWT then "--skipCaseWT " else ""
  String Pheno = if defined(PhenolyzerFile) then "--phenolyzerFile ~{PhenolyzerFile} " else ""

  String Dollar = "$"

	String totalMem = if defined(memory) then memory else memoryByThreads*threads + "M"
	Boolean inGiga = (sub(totalMem,"([0-9]+)(M|G)", "$2") == "G")
	Int memoryValue = sub(totalMem,"([0-9]+)(M|G)", "$1")
	Int totalMemMb = if inGiga then memoryValue*1024 else memoryValue
	Int memoryByThreadsMb = floor(totalMemMb/threads)

  command <<<
  set -exo pipefail

  Pheno=""
  if [[ -f "~{PhenolyzerFile}" ]]; then
    Pheno="--phenolyzerFile ~{PhenolyzerFile}"
  fi
  if [[ "~{poolSample}" != "" ]]; then
    pool="~{poolSample} \"~{sep=',' PooledSamples}\" "
  else
    pool=""
  fi

  ~{PerlExe} ~{AchabExe} \
    --vcf ~{mpavcf} \
    --outDir ~{OutDir}/ \
    --outPrefix ~{SampleID} \
    --case "~{CaseSample}" \
    ~{Dad} \
    ~{Mum} \
    ~{hope} \
    ~{candidates} \
    ~{Dollar}{Pheno} \
    --popFreqThr "~{AllelicFrequency}" \
    --filterList "~{FilterList}" \
    ~{cngGL} \
    ~{customVcf} \
    --mozaicRate "~{MozaicRate}" \
    --mozaicDP "~{MozaicDP}" \
    ~{customInfoList} \
    ~{affected} \
    ~{favGenRef} \
    ~{filtCustVcf} \
    ~{filtCustVcfReg} \
    ~{Dollar}{pool} \
    ~{sampSub} \
    ~{addCasDep} \
    ~{interVcf} \
    ~{poorCov} \
    ~{SkipCase}

 >>>

 output {
  File outAchab = "~{OutDir}/~{SampleID}_achab_catch.xlsx"
 }

 runtime {
   cpu: "~{threads}"
   requested_memory_mb_per_core: "${memoryByThreadsMb}"
 }
}
