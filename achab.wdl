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

 File AchabExe
 File OutMpa
 String? OutPhenolyzer
 String CustomVCF
 String GenesOfInterest
 Float AllelicFrequency
 Float? MozaicRate
 Float? MozaicDP
 String? CnvGeneList
 String FilterList
 String FatherSample
 String CaseSample
 String MotherSample
 String WorkflowType
 String CheckTrio
 String CustomInfo
 String SampleID
 String OutDir
 String PerlPath
 String? Affected
 String? FavouriteGeneRef
 String? FilterCustomVCF
 String? FilterCustomVCFRegex
 String? PooledSamples
 String? IDSNP
 String? SampleSubset
 String? AddCaseDepth
 File? IntersectVCF
 File? PoorCoverageFile
 File? Genemap2File
 String NewHope

 #runtime attributes
 Int Cpu
 Int Memory

 command <<<
   "${PerlPath}" "${AchabExe}" \
  --vcf "${OutMpa}" \
  --outDir "${OutDir}${SampleID}/${WorkflowType}/achab_excel/" \
  --outPrefix "${SampleID}" \
  --case "${CaseSample}" \
  --dad "${FatherSample}" \
  --mum "${MotherSample}" \
  ${CheckTrio} \
  ${NewHope} \
  --candidates "${GenesOfInterest}" \
  --phenolyzerFile "${OutPhenolyzer}" \
  --popFreqThr "${AllelicFrequency}" \
  --filterList "${FilterList}" \
  --cnvGeneList "${CnvGeneList}" \
  --customVCF "${CustomVCF}" \
  --mozaicRate "${MozaicRate}" \
  --mozaicDP "${MozaicDP}" \
  --customInfoList "${CustomInfo}" \
  
  --affected "${Affected}" \
  --favouriteGeneRef "${FavouriteGeneRef}" \
  --filterCustomVCF "${FilterCustomVCF}" \
  --filterCustomVCFRegex "${FilterCustomVCFRegex}" \
  --pooledSamples "${PooledSamples}" \
  --IDSNP "${IDSNP}" \
  --sampleSubset "${SampleSubset}" \
  ${AddCaseDepth} \
  --intersectVCF "${IntersectVCF}" \
  --poorCoverageFile "${PoorCoverageFile}" \
  --genemap2File "${Genemap2File}" 
 >>>
 
 output {
  File outAchab = "${OutDir}${SampleID}/${WorkflowType}/achab_excel/${SampleID}_achab_catch.xlsx"
 }
 runtime {                                                                                                                                                                    
  cpu: "${Cpu}"
  requested_memory_mb_per_core: "${Memory}"
 }
}
