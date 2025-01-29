# _Microdipodops_ Revision

**Citation**: Hafner, J.C., Upham, N.S., Gowen-Huang, F., and Light, J.E. Cryptic species and taxonomic revision of kangaroo mice, the rodent genus _Microdipodops_. _Journal of Mammalogy_.


## Table of contents
data
- alignmentsDNA
	- MicroAll_outPerog_allGenes.nex: full matrix in nexus format
	- MicroAll_outPerog_allGenes.phy: full matrix in phylip format

- phenotypicData
	- Hafner et al. MS cranial morphometrics 2025.xls: Raw morphometric data for 441 specimens across the original 20 variables. As explained in the text, our discriminant function analyses were based on 387 specimens (337/441 specimens had complete cranial data for the 14 variables used in our study).  
	- Hafner et al. MS colorimetrics_2025.xls: Raw colorimetric data for 678 specimens across 8 variables. As explained in the text, our discriminant function analyses were based on 678 specimens and using only 3 variables.  

finalFigures
- mainFigures
- supplementaryFigures

outputs
- phylogeneticAnalyses
	- allBayesianAnalyses_fullRuns
		- Micros_allGenes_50M_1Exp: all genes, 50 million generations, chrongram - 1 fossil calibration of exponential
		- Micros_allGenes_50M_1Uni: all genes, 50 million generations, chrongram - 1 fossil calibration of uniform
		- Micros_allGenes_50M_phylo: all genes, 50 million generations, phylogram
		- Micros_allGenes_50M_phylo_noHistorical: all genes, 50 million generations, phylogram, excluding historical specimens
		- Micros_mtDNA-only_50M_phylo: mtDNA only, 50 million generations, phylogram

	- consensusTrees: .tre and .pdf files for the MCC consensus trees of each analysis
	- plotMCC_Micros.R: R code for plotting these MCC phylogenies

- rangeMaps
	- backgroundShapes: coastlines, state borders needed to reproduce the map figures
	- RangeMap_SHPs: expert range maps for the 6 _Microdipodops_ species, and for the genus-level distribution
