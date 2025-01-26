library(ape)
library(phytools); library(phylotate)
library(phyloch) #needs to be installed from the website, not CRAN
#install.packages("phylotate")
#remotes::install_github("fmichonneau/phyloch")
#setwd("insert your directory here") 
setwd("/Users/unathan/Dropbox\ (ASU)/PROJECTS/Colllaborations/JohnJessica_Microdipodops")


#=============
# Main text FIGS
#   (1) == #4 phylo
#	(3) == #4 chrono/dated
#=============
#	(2) == #5 phylo
# LOAD PHYLOGRAM (and get the node support values + 95% HPDs per tree)
#=================================================
	treePhylo<-read_annotated(filename="final_reRuns_MB_dated-n-phylograms/infile_figtree.con_Micros_phylo_labelsX_addSpace2x.tre", format="nexus")

	posterior<-c(); height_95_HPD_MIN<-c(); height_95_HPD_MAX<-c(); height_mean<-c()
	for(j in 1:length(treePhylo$node.comment)){
		posterior[j]<-as.numeric(strsplit(strsplit(treePhylo$node.comment[j],",")[[1]][1],"=")[[1]][2])
	}
	treePhylo$posterior <- posterior[(length(treePhylo$tip.label)+1):length(treePhylo$node.comment)]

# Front matter
######
	treeToPlot<-treePhylo

# PLOT -- single-page PDFs
#=====================	

	## (1) Just M. PSAMMOPHILOMYS 
	###
	pdf(file=paste0("FIGURES/PHYLO/phyloFig1_M-Psammophilomys_spaces.pdf"), width=8.5, height=11, onefile=TRUE)

	speciesMP<-as.character(na.omit(treeToPlot$tip.label[ getDescendants(treeToPlot,getMRCA(treeToPlot,c("Mina__MLZ_1782","Alamo__MSB_35536")) ) ]))
	speciesMM<-as.character(na.omit(treeToPlot$tip.label[ getDescendants(treeToPlot,getMRCA(treeToPlot,c("Izenhood__MVZ_70918_X","Duckwater__MLZ_1997")) ) ]))

	tree<-drop.tip(treeToPlot, c("Perognathus_flavus", treeToPlot$tip.label[grep("Dipodomys",(treeToPlot$tip.label))], speciesMM) )

	#correctPP_tree3<-cbind(178:353,treeToPlot$posterior[21:196])
	correctPP_tree4<-cbind( ((178+129):353),treeToPlot$posterior[(21+129):196])
		tipColors<-rep("black",length(tree$tip.label))

	#change tips to put in parentheses
		sepNames<-do.call(rbind,strsplit(tree$tip.label,"__"))
		sepNames2<-do.call(rbind,strsplit(sepNames[,2],"_"))
		#sepNames[which(sepNames[,5]!="X"),5] <- ""
		#sepNames[which(sepNames[,5]=="X"),5] <- "*"
		newTips<-paste0(sepNames[,1],"__(",sepNames2[,1],"_",sepNames2[,2],")")
	
	tree2<-tree
	tree2$tip.label<-newTips

	#quartz(width=8.5, height=22)
	plot(ladderize(tree2), cex=0.6, label.offset=0.00005, tip.color=tipColors, y.lim=c(-15,50), font=1) #x.lim=c(-7,40))#
	#node.support(tree$posterior, mode="dots", col = "red", cex=0.3)
	#node.support(tree$posterior, mode="dots", cutoff=0.95, col = "black", cex=0.55)
	#node.support(tree$posterior, mode="numbers", digits=2,pos="below", col = "black", font=2, cex=0.5)
	#node.support(correctPP_tree4[,2], mode="dots", col = "red", cex=0.3)
	node.support(correctPP_tree4[,2], mode="dots", cutoff=0.90, col = "black", cex=0.55)
	node.support(correctPP_tree4[,2], mode="numbers", cutoff=0.90, digits=2,pos="below", col = "black", font=2, cex=0.5)

	add.scale.bar(x=par()$usr[1],y=1,lwd=1, cex=0.7)
	text(x=0.003,y=1,labels="expected substitutions per site", adj=0, cex=0.7)

	mtext(side=3, text="phylogeny Fig #1", font=2, line= 0, adj=0, cex=0.7)
	mtext(side=3, text="M. (Psammophilomys)", font=3, line= -1, adj=0, cex=0.7)

	dev.off()

	## (2) Just M. MICRODIPODOPS 
	###
	pdf(file=paste0("FIGURES/PHYLO/phyloFig2_M-Microdipodops_spaces.pdf"), width=8.5, height=11, onefile=TRUE)

	speciesMP<-as.character(na.omit(treeToPlot$tip.label[ getDescendants(treeToPlot,getMRCA(treeToPlot,c("Mina__MLZ_1782","Alamo__MSB_35536")) ) ]))
	speciesMM<-as.character(na.omit(treeToPlot$tip.label[ getDescendants(treeToPlot,getMRCA(treeToPlot,c("Izenhood__MVZ_70918_X","Duckwater__MLZ_1997")) ) ]))

	tree<-drop.tip(treeToPlot, c("Perognathus_flavus", treeToPlot$tip.label[grep("Dipodomys",(treeToPlot$tip.label))], speciesMP) )

	#correctPP_tree3<-cbind(178:353,treeToPlot$posterior[21:196])
	#correctPP_tree4<-cbind( ((178+129):353),treeToPlot$posterior[(21+129):196])
	correctPP_tree5<-treeToPlot$posterior[22:(21+128)]
	#tree$posterior<-correctPP_tree5[,2]
		tipColors<-rep("black",length(tree$tip.label))

	#change tips to put in parentheses
		sepNames<-do.call(rbind,strsplit(tree$tip.label,"__"))
		sepNames2<-do.call(rbind,strsplit(sepNames[,2],"_"))
		#sepNames[which(sepNames[,5]!="X"),5] <- ""
		#sepNames[which(sepNames[,5]=="X"),5] <- "*"
		newTips<-paste0(sepNames[,1],"__(",sepNames2[,1],"_",sepNames2[,2],")")
	
	tree2<-tree
	tree2$tip.label<-newTips

	par(xpd=NA) #l value or ‘NA’.  If ‘FALSE’, all plotting is
     #     clipped to the plot region, if ‘TRUE’, all plotting is
      #    clipped to the figure region, and if ‘NA’, all plotting is
       #   clipped to the device region.  See also ‘clip’.

	#quartz(width=8.5, height=22)
	plot(ladderize(tree2), cex=0.47, label.offset=0.00005, tip.color=tipColors, y.lim=c(14,125), font=1) #x.lim=c(-7,40))#
	#node.support(tree$posterior, mode="dots", col = "red", cex=0.3)
	#node.support(tree$posterior, mode="dots", cutoff=0.95, col = "black", cex=0.55)
	#node.support(tree$posterior, mode="numbers", digits=2,pos="below", col = "black", font=2, cex=0.4)
	#node.support(correctPP_tree5, mode="dots", col = "red", cex=0.3)
	node.support(correctPP_tree5, mode="dots", cutoff=0.90, col = "black", cex=0.55)
	node.support(correctPP_tree5, mode="numbers", cutoff=0.90, digits=2,pos="below", col = "black", font=2, cex=0.5)

	add.scale.bar(x=par()$usr[1],y=1,lwd=1, cex=0.7)
	text(x=0.007,y=1,labels="expected substitutions per site", adj=0, cex=0.7)

	mtext(side=3, text="phylogeny Fig #2", font=2, line= 0, adj=0, cex=0.7)
	mtext(side=3, text="M. (Microdipodops)", font=3, line= -1, adj=0, cex=0.7)

	dev.off()

	## (1,2 SUB) the inset Fig of whole genus
	###
	tree<-drop.tip(treeToPlot, c("Perognathus_flavus", treeToPlot$tip.label[grep("Dipodomys",(treeToPlot$tip.label))]) )
	correctPP_tree3<-cbind(178:353,treeToPlot$posterior[21:196])

	pdf(file=paste0("FIGURES/PHYLO/phyloFig12sub_genusSmall.pdf"), width=6, height=8, onefile=TRUE)

	#quartz(width=8.5, height=22)
	plot(ladderize(tree), cex=0.35, label.offset=0.00005, show.tip.label=FALSE, y.lim=c(5,180)) #x.lim=c(-7,40))#

	dev.off()


# LOAD CHRONOGRAM (and get the node support values + 95% HPDs per tree)
#=================================================
	treeExp<-read_annotated(filename="final_reRuns_MB_dated-n-phylograms/infile_figtree.con_Micros_1Exp_labelsX_addSpace.tre", format="nexus")
	posterior<-c(); height_95_HPD_MIN<-c(); height_95_HPD_MAX<-c(); height_mean<-c()
	for(j in 1:length(treeExp$node.comment)){
		posterior[j]<-as.numeric(strsplit(strsplit(treeExp$node.comment[j],",")[[1]][1],"=")[[1]][2])
		height_95_HPD_MIN[j]<-as.numeric(strsplit(strsplit(treeExp$node.comment[j],",")[[1]][9],"=")[[1]][2])
		height_95_HPD_MAX[j]<-as.numeric(strsplit(treeExp$node.comment[j],",")[[1]][10])
		height_mean[j]<-as.numeric(strsplit(strsplit(treeExp$node.comment[j],",")[[1]][7],"=")[[1]][2])
	}
	treeExp$posterior <- posterior[(length(treeExp$tip.label)+1):length(treeExp$node.comment)]
	treeExp$'height_95%_HPD_MIN' <- height_95_HPD_MIN[(length(treeExp$tip.label)+1):length(treeExp$node.comment)]
	treeExp$'height_95%_HPD_MAX' <- height_95_HPD_MAX[(length(treeExp$tip.label)+1):length(treeExp$node.comment)]
	treeExp$height_mean <- height_mean[(length(treeExp$tip.label)+1):length(treeExp$node.comment)]

# Front matter
######
	treeToPlot<- treeExp

	## (3) CHRONOGRAM -- Major clades + DIPODOMYS 
	###
	# keep only 1 rep per clade
	toKeep<-c(treeExp$tip.label[grep("Dipodomys",(treeExp$tip.label))], "Mina__MLZ_1782","Alamo__MSB_35536", "Owyhee__MLZ_2180", "FortRock__MLZ_2174","Minersville__MLZ_2071","WEureka__MLZ_2031")
	toDrop<-setdiff(treeExp$tip.label,toKeep)
	treeExp_simp <- drop.tip2(treeExp,toDrop)
		# relabel:
		oldLabel<-c("Mina__MLZ_1782","Alamo__MSB_35536", "Owyhee__MLZ_2180", "FortRock__MLZ_2174","Minersville__MLZ_2071","WEureka__MLZ_2031")
		newLabel<-c("M. (P.) pallidus", "M. (P.) ruficollaris", "M. (M.) megacephalus", "M. (M.) oregonus", "M. (M.) albiventer", "M. (M.) polionotus")
		for(i in 1: length(oldLabel)){
			tipNum<-match(oldLabel[i], treeExp_simp$tip.label)
			treeExp_simp$tip.label[tipNum]<-newLabel[i]
		}

	# PLOT simplified version...
	pdf(file=paste0("FIGURES/PHYLO/phyloFig3_Micros-vs-Dipos.pdf"), width=8.5, height=11, onefile=TRUE)

	tree<-treeExp_simp

	x1= -7 ; x2= 30
	y1= -15 ; y2= 30

	# plot dummy tree
	plot(ladderize(tree), cex=0.8, label.offset=0.4, tip.color="white", edge.width=2, y.lim=c(y1,y2), x.lim=c(x1,x2))
		# put rectangle
		root<-max(branching.times(tree))
		rect(xleft=root-2.6, ybottom=0, xright=root-0.8, ytop=27, col=grey(0.5, alpha=0.3), border=NA)

	par(new=TRUE)
	#quartz(width=8.5, height=22)
	plot(ladderize(tree), cex=0.8, label.offset=0.4, tip.color="black", edge.width=2, y.lim=c(y1,y2), x.lim=c(x1,x2))

	HPDbars(tree, label="height_95%_HPD", broken=T, lwd=3, col=hsv(0.65,1,1,alpha=0.7))
	node.support(tree$posterior, mode="dots", col = "red", cex=0.4)
	node.support(tree$posterior, mode="dots", cutoff=0.90, col = "black", cex=0.75)
	node.support(branching.times(tree), mode="numbers", digits=1,pos="above", col = "black", font=2, cex=0.8)
	#node.support(tree$posterior, mode="numbers", digits=2,pos="below", col = "grey50", font=2, cex=0.8)

	data(gradstein04)
	axisGeo(GTS = gradstein04, unit = c("epoch","period"), cex = 0.8, ages=FALSE, gridty=3, gridcol="grey50")
	axisPhylo(cex.axis=0.8, pos=-2)#, mgp=c(0,0.1,0))

	mtext(side=3, text="phylogeny Fig #3", font=2, line= 1, adj=0, cex=0.7)
	mtext(side=3, text="Node support: black if PP >= 0.90; red if PP < 0.90", font=1, line= 0, adj=0, cex=0.7)
	mtext(side=3, text="6 lineages of Microdipodops vs. 20 species of Dipodomys (gray box is 95% HPD of Microdipodops sister species)", font=1, line= -1, adj=0, cex=0.7)

	dev.off()


#=============
# Supplemental FIGS
#   (6) == GENE TREES 
#=============
# LOAD EACH TREE -- RAxML -- mtDNA, GHR, IRBP, EN2 -- (and get the node support values + 95% HPDs per tree)
#=================================================
	GENES<-c("MtDNA", "GHR", "IRBP", "EN2")

	treePhylo_GENES<-list()
	for(i in 1:length(GENES)){
		treePhylo<-read_annotated(filename=paste0("geneTrees/Micro",GENES[i],"_11Feb19_RAxML_1000bs.nex"), format="nexus")

		posterior<-c(); height_95_HPD_MIN<-c(); height_95_HPD_MAX<-c(); height_mean<-c()
		for(j in 1:length(treePhylo$node.comment)){
			#posterior[j]<-as.numeric(strsplit(strsplit(treePhylo$node.comment[j],"=")[[1]][1],"=")[[1]][2])
			posterior[j]<-as.numeric( strsplit(treePhylo$node.comment[j],"=")[[1]][2] )
		}
		treePhylo$posterior <- posterior[(length(treePhylo$tip.label)+1):length(treePhylo$node.comment)]
	
	treePhylo_GENES[[i]]<-treePhylo
	}

	## get the node membership for ALL NODES
	##nodesGrThan50 <- treePhylo$posterior[which(treePhylo$posterior>50)]
	#descendantsPerNode<-list()
	#for(q in 1:length(treePhylo$posterior)){
	#	node<-length(treePhylo$posterior)+q
	#	descendantsPerNode[[q]]<-as.character(na.omit(treePhylo$tip.label[ getDescendants(treePhylo, node) ]))
	#
	#}
	#new<-do.call(rbind,descendantsPerNode)

# get all the taxa, assign them to clades:
	allTips<-unique(c(treePhylo_GENES[[1]]$tip.label, treePhylo_GENES[[2]]$tip.label, 
						treePhylo_GENES[[3]]$tip.label, treePhylo_GENES[[4]]$tip.label)) # 178 total

# Load in the full clade members (6 species)
	cladeMembership<-read.csv(file="MicrodipodopsSampling12March2021_v3_2023_geneTrees.csv")
		# 178 tips now
	#cladeMembership<-read.csv(file="MicrodipodopsSampling12March2021_v2_2023.csv")
		# 176 species originally...
		#	# which were not present?
		#	dd<-left_join(x=cbind.data.frame(allTips=tree$tip.label,rep(0,length(allTips))), y=cladeMembership, by=c("allTips" = "TaxonName"))
		#	dd[is.na(dd[,4]),]
		#	# SOLVED THESE
		#		#                  allTips rep(0, length(allTips)) Clade Species
		#		# 147       GerlachMLZ2097                       0  <NA>    <NA>
		#		# 153       GerlachMLZ2089                       0  <NA>    <NA>
		#		# 155       GerlachMLZ2092                       0  <NA>    <NA>
		#		# 157       GerlachMLZ2099                       0  <NA>    <NA>
		#		# 158 NWGolcondaMVZ70842_X                       0  <NA>    <NA>
		#		# 159       GerlachMLZ2094                       0  <NA>    <NA>
		#		# 178     Dipodomys_agilis                       0  <NA>    <NA>
	# load in number translation table
		# localityList<-read.csv(file="MicrodipodopsSampling12March2021_v3_localityList.csv")
		# joined<-left_join(cladeMembership, localityList, by= c("LocalityName" = "General.Locality.noSpaces"))
		# write.csv(joined[,1:9], file="MicrodipodopsSampling12March2021_v3_2023_updatedNums.csv")

# Front matter
######

	library(viridis)
	# PLOT -- all GENE TREES together (separate pages)
	#=====================	
	YLIMS<-list(c(5,180), c(5,140), c(5,140), c(5,140))
	num1<-0
	XLIMS<-list(c(num1,0.5), c(num1,0.05), c(num1,0.06), c(num1,0.09))
	
	SCALES<-c(0.08,0.008,0.008,0.015)
	OFFSET<-c(0.004,0.0004,0.0004,0.0012)
#	colors<-c(magma(12)[c(3,5,8,10)],viridis(12)[c(5,10)])
#	colors<-c(magma(12)[c(2,4,7,9)],viridis(5)[4],"goldenrod1")#viridis(12)[c(1,3,5,7,9,11)]
	colors<-c("#d73027","#231F20","#fc8d59","#52bfdbff","#939598","#542788")#"#7f5da6bf")
#	colors<- rainbow(6) #make this a rainbow...
	PART<-c("(A)","(B)","(C)","(D)" )

#	pdf(file=paste0("FIGURES/PHYLO/geneTrees_ALL_cut75_colorBranches_wLocNumsCorrect_narrowLong_contrastReady.pdf"), width=8.5, height=20, onefile=TRUE)
	pdf(file=paste0("FIGURES/PHYLO/SD4_geneTrees_ALL_cut75_colorBranches_wLocNames_narrowLong_Nov2024_v4.pdf"), width=8.5, height=20, onefile=TRUE)

	for(i in 1:length(treePhylo_GENES)){
		tree<-treePhylo_GENES[[i]]
			tipColors<-rep("black",length(tree$tip.label))

			# set tip colors by clade
			for(j in 1:length(tree$tip.label)){
				tip <- tree$tip.label[j]
				species <- cladeMembership[which(cladeMembership$TaxonName==tip), "Species"]
				if(species == "albiventer"){
					tipColors[j] <- colors[1]
				} else if(species == "megacephalus"){
					tipColors[j] <- colors[2]
				} else if(species == "oregonus"){
					tipColors[j] <- colors[3]
				} else if(species == "polionotus"){
					tipColors[j] <- colors[4]
				} else if(species == "pallidus"){
					tipColors[j] <- colors[5]
				} else if(species == "ruficollaris"){
					tipColors[j] <- colors[6]
				} else { tipColors[j] <- "black" }
			}

		# change locality names to NUMBERS
	#		treeReorder<-match(tree$tip.label, cladeMembership$TaxonName)
	#		dat_reordered<-cladeMembership[treeReorder,]
	#		newTips<-paste0(dat_reordered$LocalityNum,": ",dat_reordered$MuseumNum)

		# change locality names to BETTER FORMAT
			treeReorder<-match(tree$tip.label, cladeMembership$TaxonName)
			dat_reordered<-cladeMembership[treeReorder,]
			newTips<-paste0(dat_reordered$LocalityName,": ",dat_reordered$MuseumNum)

		# change aDNA tips
	#	 	sepNames<-do.call(rbind,strsplit(tree$tip.label,"_"))
	#		sepNames[which(sepNames[,5]!="X"),5] <- ""
	#		sepNames[which(sepNames[,5]=="X"),5] <- "*"
	#		newTips<-paste0(sepNames[,1],"_",sepNames[,2],"_(",sepNames[,3],"_",sepNames[,4],")_",sepNames[,5])
		
		tree2<-tree
		tree2$tip.label<-newTips

		plot(ladderize(tree2), type="phylogram",cex=0.7, label.offset=OFFSET[i], tip.color=tipColors, font=2, y.lim=YLIMS[[i]], x.lim=XLIMS[[i]])#
	#	node.support(tree$posterior, mode="dots", col = "red", cex=0.3)
	#	node.support(tree$posterior, mode="dots", cutoff=50, col = gray(0.5), cex=0.55)
	#	node.support(tree$posterior, mode="numbers", cutoff=50, digits=2,pos="below", col = gray(0.5), font=2, cex=0.6)
		node.support(tree$posterior, mode="dots", cutoff=75, col = "black", cex=1.3)
		node.support(tree$posterior, mode="numbers", cutoff=75, digits=2,pos="below", col = "black", font=2, cex=0.85)

		add.scale.bar(x=par()$usr[1],y=1,lwd=1, cex=0.7)
		text(x=SCALES[i],y=1,labels="expected substitutions per site", adj=0, cex=0.7)

		mtext(side=3, text=paste0(PART[i]), font=2, line= 1, adj=0, cex=0.7)
		mtext(side=3, text=paste0("gene tree ", GENES[i]), font=2, line= 0, adj=0, cex=0.7)
		mtext(side=3, text="Microdipodops", font=3, line= -1, adj=0, cex=0.7)

		legend(x=par()$usr[1],y=50,legend=c("M.(M.) albiventer",  "M.(M.) megacephalus", "M.(M.) oregonus", "M.(M.) polionotus",  "M.(P.) pallidus", "M.(P.) ruficollaris"), 
			col=colors, pch=15, text.font=3, pt.cex=2)

	}
	dev.off()



# PHYLOGRAMS:
#================

#=============
# Supplemental FIGS (phylograms)
#   == #1, #2, #3, #4, #5 
#=============
# LOAD TREE (and get the node support values + 95% HPDs per tree)
#=================================================

	#treePhylo<-read_annotated(filename="final_reRuns_MB_dated-n-phylograms/infile_figtree.con_Micros_phylo_labelsX_addSpace.tre", format="nexus")
#	treePhylo<-read_annotated(filename="final_reRuns_MB_dated-n-phylograms/infile_figtree.con_Micros_phylo_mtDNA-only.tre", format="nexus")
	treePhylo<-read_annotated(filename="final_reRuns_MB_dated-n-phylograms/infile_figtree.con_Micros_phylo_noHistorical.tre", format="nexus")

	posterior<-c(); height_95_HPD_MIN<-c(); height_95_HPD_MAX<-c(); height_mean<-c()
	for(j in 1:length(treePhylo$node.comment)){
		posterior[j]<-as.numeric(strsplit(strsplit(treePhylo$node.comment[j],",")[[1]][1],"=")[[1]][2])
	}
	treePhylo$posterior <- posterior[(length(treePhylo$tip.label)+1):length(treePhylo$node.comment)]

# Front matter
######
#	treeType<-"Phylogram (no calibrations): MrBayes Allcompat consensus of 10k trees"
#	treeType<-"** mtDNA only ** Phylogram (no calibrations): MrBayes Allcompat consensus of 10k trees"
	treeType<-"** all genes; excluding historical ** Phylogram (no calibrations): MrBayes Allcompat consensus of 10k trees"
	cal<-"phylogram-Allcompat"
	treeToPlot<-treePhylo

# PLOT -- multi-page PDF with each of the iterations
#=====================	

#	pdf(file=paste0("FIGURES/PHYLO/MCCtree_",cal,"_plot1-5_2022_ALL_v2.pdf"), width=8.5, height=11, onefile=TRUE)
#	pdf(file=paste0("FIGURES/PHYLO/SD5_MCCtree_",cal,"_plot1_allTaxa_v4_mtDNAonly.pdf"), width=8.5, height=11, onefile=TRUE)
	pdf(file=paste0("FIGURES/PHYLO/SD6_MCCtree_",cal,"_plot1_allGenes_v4_noHist.pdf"), width=8.5, height=11, onefile=TRUE)

	## (1) ALL SPECIES
	###
	#pdf(file=paste0("FIGURES/PHYLO/MCCtree_plot_Microdipodops_phylogram_noAges_PPdots_ColxTips_2022_ALL.pdf"), width=20, height=40, onefile=TRUE)
	#pdf(file=paste0("FIGURES/PHYLO/MCCtree_",cal,"_plot1_2022_allTips_v2.pdf"), width=8.5, height=11, onefile=TRUE)

	tree<-treeToPlot
	#get tip colors
		sepNames<-do.call(rbind,strsplit(tree$tip.label,"_"))
		tipColors<-rep("black",length(tree$tip.label))
	#	tipColors[which(sepNames[,5]=="X")] <- "red"

	#quartz(width=8.5, height=22)
#	plot(ladderize(tree), cex=0.3, label.offset=0.001, tip.color=tipColors, y.lim=c(5,200))#, x.lim=c(-7,40))#, y.lim=c(130,3950))
	plot(ladderize(tree), cex=0.3, label.offset=0.001, tip.color=tipColors, y.lim=c(0,165))#, x.lim=c(-7,40))#, y.lim=c(130,3950))
	node.support(tree$posterior, mode="dots", col = "red", cex=0.25)
	node.support(tree$posterior, mode="dots", cutoff=0.95, col = "black", cex=0.55)
	node.support(tree$posterior, mode="numbers", digits=2,pos="below", col = "grey50", font=2, cex=0.25)

	add.scale.bar(x=par()$usr[1],y=-1.5,lwd=1, cex=0.7)
#	text(x=0.065,y=1,labels="expected substitutions per site", adj=0, cex=0.7)
	text(x=0.065,y=-1.5,labels="expected substitutions per site", adj=0, cex=0.7)

	mtext(side=3, text=treeType, font=2, line=0, adj=0, cex=0.7)
#	mtext(side=3, text="Plot #1 - Node support: black if PP >= 0.95; red if PP < 0.95 - Tips: red if ancient DNA", font=2, line= -1, adj=0, cex=0.7)
#	mtext(side=3, text="Node support: black if PP >= 0.95; red if PP < 0.95 - Tips: red if ancient DNA", font=2, line= -1, adj=0, cex=0.7)
	mtext(side=3, text="Node support: black if PP >= 0.95; red if PP < 0.95", font=2, line= -1, adj=0, cex=0.7)
	#mtext(side=3, text="197 total tips (176 Microdipodops), 20 Dipodomys, 1 Perognathus outgroup", font=3, line= -2, adj=0, cex=0.7)
	mtext(side=3, text="all Microdipodops, 20 Dipodomys, 1 Perognathus outgroup", font=3, line= -2, adj=0, cex=0.7)

dev.off()


	## (2) MINUS OUTGROUP 
	###
	#pdf(file=paste0("FIGURES/PHYLO/MCCtree_",cal,"_plot2_2022_noOut_v2.pdf"), width=8.5, height=11, onefile=TRUE)

	#tree<-treeExp
	tree<-drop.tip(treeToPlot,"Perognathus_flavus")
	#get tip colors
		sepNames<-do.call(rbind,strsplit(tree$tip.label,"_"))
		tipColors<-rep("black",length(tree$tip.label))
		tipColors[which(sepNames[,5]=="X")] <- "red"

	#quartz(width=8.5, height=22)
	plot(ladderize(tree), cex=0.35, label.offset=0.0005, tip.color=tipColors, y.lim=c(5,200))#, x.lim=c(-7,40))#, y.lim=c(130,3950))
	node.support(tree$posterior, mode="dots", col = "red", cex=0.3)
	node.support(tree$posterior, mode="dots", cutoff=0.95, col = "black", cex=0.55)
	node.support(tree$posterior, mode="numbers", digits=2,pos="below", col = "grey50", font=2, cex=0.3)

	add.scale.bar(x=par()$usr[1],y=1,lwd=1, cex=0.7)
	text(x=0.025,y=1,labels="expected substitutions per site", adj=0, cex=0.7)

	mtext(side=3, text=treeType, font=2, line=0, adj=0, cex=0.7)
	mtext(side=3, text="Plot #2 - Node support: black if PP >= 0.95; red if PP < 0.95 - Tips: red if ancient DNA", font=2, line= -1, adj=0, cex=0.7)
	#mtext(side=3, text="196 total tips (176 Microdipodops), 20 Dipodomys, excluding outgroup", font=3, line= -2, adj=0, cex=0.7)
	mtext(side=3, text="All Microdipodops, 20 Dipodomys, excluding outgroup", font=3, line= -2, adj=0, cex=0.7)

	#dev.off()


	## (3) MINUS DIPODOMYS 
	###
	#pdf(file=paste0("FIGURES/PHYLO/MCCtree_",cal,"_plot3_2022_noOutNoDipos_v2.pdf"), width=8.5, height=11, onefile=TRUE)

	tree<-drop.tip(treeToPlot, c("Perognathus_flavus", treeToPlot$tip.label[grep("Dipodomys",(treeToPlot$tip.label))]) )
	#correctPP_fullTree<-cbind(199:394,treeToPlot$posterior)
	#correctPP_tree2<-cbind(198:393,tree2$posterior)
	correctPP_tree3<-cbind(178:353,treeToPlot$posterior[21:196])

	#get tip colors
		sepNames<-do.call(rbind,strsplit(tree$tip.label,"_"))
		tipColors<-rep("black",length(tree$tip.label))
		tipColors[which(sepNames[,5]=="X")] <- "red"

	#quartz(width=8.5, height=22)
	plot(ladderize(tree), cex=0.35, label.offset=0.00005, tip.color=tipColors, y.lim=c(5,180)) #x.lim=c(-7,40))#
	#node.support(tree$posterior, mode="dots", col = "red", cex=0.5)
	#node.support(tree$posterior, mode="dots", cutoff=0.95, col = "black", cex=0.55)
	#node.support(tree$posterior, mode="numbers", digits=2,pos="below", col = "black", font=2, cex=0.7)
	node.support(correctPP_tree3[,2], mode="dots", col = "red", cex=0.3)
	node.support(correctPP_tree3[,2], mode="dots", cutoff=0.95, col = "black", cex=0.55)
	node.support(correctPP_tree3[,2], mode="numbers", digits=2,pos="below", col = "grey50", font=2, cex=0.3)

	add.scale.bar(x=par()$usr[1],y=1,lwd=1, cex=0.7)
	text(x=0.0075,y=1,labels="expected substitutions per site", adj=0, cex=0.7)

	mtext(side=3, text=treeType, font=2, line=0, adj=0, cex=0.7)
	mtext(side=3, text="Plot #3 - Node support: black if PP >= 0.95; red if PP < 0.95 - Tips: red if ancient DNA", font=2, line= -1, adj=0, cex=0.7)
	#mtext(side=3, text="176 Microdipodops only, excluding Dipodomys and outgroup", font=3, line= -2, adj=0, cex=0.7)
	mtext(side=3, text="All Microdipodops only, excluding Dipodomys and outgroup", font=3, line= -2, adj=0, cex=0.7)

	#dev.off()

	## (4) Just M. PSAMMOPHILOMYS 
	###
	#pdf(file=paste0("FIGURES/PHYLO/MCCtree_",cal,"_plot4_2022_onlyMP_v2.pdf"), width=8.5, height=11, onefile=TRUE)

	speciesMP<-as.character(na.omit(treeToPlot$tip.label[ getDescendants(treeToPlot,getMRCA(treeToPlot,c("Mina__MLZ_1782","Alamo__MSB_35536")) ) ]))
	speciesMM<-as.character(na.omit(treeToPlot$tip.label[ getDescendants(treeToPlot,getMRCA(treeToPlot,c("Izenhood__MVZ_70918_X","Duckwater__MLZ_1997")) ) ]))

	tree<-drop.tip(treeToPlot, c("Perognathus_flavus", treeToPlot$tip.label[grep("Dipodomys",(treeToPlot$tip.label))], speciesMM) )

	#correctPP_tree3<-cbind(178:353,treeToPlot$posterior[21:196])
	correctPP_tree4<-cbind( ((178+129):353),treeToPlot$posterior[(21+129):196])

	#get tip colors
		sepNames<-do.call(rbind,strsplit(tree$tip.label,"_"))
		tipColors<-rep("black",length(tree$tip.label))
		tipColors[which(sepNames[,5]=="X")] <- "red"

	#quartz(width=8.5, height=22)
	plot(ladderize(tree), cex=0.6, label.offset=0.00005, tip.color=tipColors, y.lim=c(-15,50)) #x.lim=c(-7,40))#
	#node.support(tree$posterior, mode="dots", col = "red", cex=0.3)
	#node.support(tree$posterior, mode="dots", cutoff=0.95, col = "black", cex=0.55)
	#node.support(tree$posterior, mode="numbers", digits=2,pos="below", col = "black", font=2, cex=0.5)
	node.support(correctPP_tree4[,2], mode="dots", col = "red", cex=0.3)
	node.support(correctPP_tree4[,2], mode="dots", cutoff=0.95, col = "black", cex=0.55)
	node.support(correctPP_tree4[,2], mode="numbers", digits=2,pos="below", col = "grey50", font=2, cex=0.5)

	add.scale.bar(x=par()$usr[1],y=1,lwd=1, cex=0.7)
	text(x=0.003,y=1,labels="expected substitutions per site", adj=0, cex=0.7)

	mtext(side=3, text=treeType, font=2, line=0, adj=0, cex=0.7)
	mtext(side=3, text="Plot #4 - Node support: black if PP >= 0.95; red if PP < 0.95 - Tips: red if ancient DNA", font=2, line= -1, adj=0, cex=0.7)
	#mtext(side=3, text="47 M.(P) only, excluding M(M.), Dipodomys, and outgroup", font=3, line= -2, adj=0, cex=0.7)
	mtext(side=3, text="M.(P) only, excluding M(M.), Dipodomys, and outgroup", font=3, line= -2, adj=0, cex=0.7)

	#dev.off()

	## (5) Just M. MICRODIPODOPS 
	###
	#pdf(file=paste0("FIGURES/PHYLO/MCCtree_",cal,"_plot5_2022_onlyMM_v2.pdf"), width=8.5, height=11, onefile=TRUE)

	speciesMP<-as.character(na.omit(treeToPlot$tip.label[ getDescendants(treeToPlot,getMRCA(treeToPlot,c("Mina__MLZ_1782","Alamo__MSB_35536")) ) ]))
	speciesMM<-as.character(na.omit(treeToPlot$tip.label[ getDescendants(treeToPlot,getMRCA(treeToPlot,c("Izenhood__MVZ_70918_X","Duckwater__MLZ_1997")) ) ]))

	tree<-drop.tip(treeToPlot, c("Perognathus_flavus", treeToPlot$tip.label[grep("Dipodomys",(treeToPlot$tip.label))], speciesMP) )

	#correctPP_tree3<-cbind(178:353,treeToPlot$posterior[21:196])
	#correctPP_tree4<-cbind( ((178+129):353),treeToPlot$posterior[(21+129):196])
	correctPP_tree5<-treeToPlot$posterior[22:(21+128)]
	#tree$posterior<-correctPP_tree5[,2]

	#get tip colors
		sepNames<-do.call(rbind,strsplit(tree$tip.label,"_"))
		tipColors<-rep("black",length(tree$tip.label))
		tipColors[which(sepNames[,5]=="X")] <- "red"

	#quartz(width=8.5, height=22)
	plot(ladderize(tree), cex=0.47, label.offset=0.00005, tip.color=tipColors, y.lim=c(5,130)) #x.lim=c(-7,40))#
	#node.support(tree$posterior, mode="dots", col = "red", cex=0.3)
	#node.support(tree$posterior, mode="dots", cutoff=0.95, col = "black", cex=0.55)
	#node.support(tree$posterior, mode="numbers", digits=2,pos="below", col = "black", font=2, cex=0.4)
	node.support(correctPP_tree5, mode="dots", col = "red", cex=0.3)
	node.support(correctPP_tree5, mode="dots", cutoff=0.95, col = "black", cex=0.55)
	node.support(correctPP_tree5, mode="numbers", digits=2,pos="below", col = "grey50", font=2, cex=0.5)

	add.scale.bar(x=par()$usr[1],y=1,lwd=1, cex=0.7)
	text(x=0.007,y=1,labels="expected substitutions per site", adj=0, cex=0.7)

	mtext(side=3, text=treeType, font=2, line=0, adj=0, cex=0.7)
	mtext(side=3, text="Plot #5 - Node support: black if PP >= 0.95; red if PP < 0.95 - Tips: red if ancient DNA", font=2, line= -1, adj=0, cex=0.7)
	#mtext(side=3, text="129 M.(M) only, excluding M(P.), Dipodomys, and outgroup", font=3, line= -2, adj=0, cex=0.7)
	mtext(side=3, text="M.(M) only, excluding M(P.), Dipodomys, and outgroup", font=3, line= -2, adj=0, cex=0.7)

	dev.off()


# Subset PHYLOGRAMS: Excluding historical specimens, MtDNA only
#================

	#treePhylo<-read_annotated(filename="final_reRuns_MB_dated-n-phylograms/infile_figtree.con_Micros_phylo_noHistorical.tre", format="nexus")
	treePhylo<-read_annotated(filename="final_reRuns_MB_dated-n-phylograms/infile_figtree.con_Micros_phylo_mtDNA-only.tre", format="nexus")

	posterior<-c(); height_95_HPD_MIN<-c(); height_95_HPD_MAX<-c(); height_mean<-c()
	for(j in 1:length(treePhylo$node.comment)){
		posterior[j]<-as.numeric(strsplit(strsplit(treePhylo$node.comment[j],",")[[1]][1],"=")[[1]][2])
	}
	treePhylo$posterior <- posterior[(length(treePhylo$tip.label)+1):length(treePhylo$node.comment)]

# Load in the full clade members (6 species)
	cladeMembership<-read.csv(file="MicrodipodopsSampling12March2021_v3_2023_geneTrees.csv")
		# 178 tips now

		# Get list of aDNA tips to highlight red
			treeReorder<-match(treePhylo$tip.label, cladeMembership$TaxonName_noX)
			dat_reordered<-cladeMembership[treeReorder,]
			aDNA_tips<-dat_reordered[which(dat_reordered$Ancient==1),"TaxonName_noX"]

# Front matter
######
	treeType<-"** mtDNA only ** Phylogram (no calibrations): MrBayes Allcompat consensus of 10k trees"
#	treeType<-"** all genes; excluding historical ** Phylogram (no calibrations): MrBayes Allcompat consensus of 10k trees"
	cal<-"phylogram-Allcompat"
	treeToPlot<-treePhylo

# PLOT -- multi-page PDF with each of the iterations
#=====================	

	pdf(file=paste0("FIGURES/PHYLO/SD5_MCCtree_",cal,"_plot1_allTaxa_v4_mtDNAonly-redTips.pdf"), width=8.5, height=11, onefile=TRUE)

	## (1) ALL SPECIES
	###
	#pdf(file=paste0("FIGURES/PHYLO/MCCtree_plot_Microdipodops_phylogram_noAges_PPdots_ColxTips_2022_ALL.pdf"), width=20, height=40, onefile=TRUE)
	#pdf(file=paste0("FIGURES/PHYLO/MCCtree_",cal,"_plot1_2022_allTips_v2.pdf"), width=8.5, height=11, onefile=TRUE)

	tree<-treeToPlot
	#get tip colors
		tipColors<-rep("black",length(tree$tip.label))
		tipColors[match(aDNA_tips, tree$tip.label)] <- "red"
		#datz<-cbind.data.frame(tree$tip.label,tipColors)

	#quartz(width=8.5, height=22)
	plot(ladderize(tree), cex=0.3, label.offset=0.001, tip.color=tipColors, y.lim=c(5,200))#, x.lim=c(-7,40))#, y.lim=c(130,3950))
	node.support(tree$posterior, mode="dots", col = "red", cex=0.25)
	node.support(tree$posterior, mode="dots", cutoff=0.95, col = "black", cex=0.55)
	node.support(tree$posterior, mode="numbers", digits=2,pos="below", col = "grey50", font=2, cex=0.25)

	add.scale.bar(x=par()$usr[1],y=1,lwd=1, cex=0.7)
	VAR=0.065	
	text(x=VAR,y=1,labels="expected substitutions per site", adj=0, cex=0.7)

	mtext(side=3, text=treeType, font=2, line=0, adj=0, cex=0.7)
	mtext(side=3, text="Node support: black if PP >= 0.95; red if PP < 0.95 - Tips: red if historical DNA samples", font=2, line= -1, adj=0, cex=0.7)
	#mtext(side=3, text="197 total tips (176 Microdipodops), 20 Dipodomys, 1 Perognathus outgroup", font=3, line= -2, adj=0, cex=0.7)
	mtext(side=3, text="all Microdipodops, 20 Dipodomys, 1 Perognathus outgroup", font=3, line= -2, adj=0, cex=0.7)

dev.off()


	## (2) MINUS OUTGROUP 
	###
	#pdf(file=paste0("FIGURES/PHYLO/MCCtree_",cal,"_plot2_2022_noOut_v2.pdf"), width=8.5, height=11, onefile=TRUE)

	#tree<-treeExp
	tree<-drop.tip(treeToPlot,"Perognathus_flavus")
	#get tip colors
		tipColors<-rep("black",length(tree$tip.label))
		tipColors[match(aDNA_tips, tree$tip.label)] <- "red"

	#quartz(width=8.5, height=22)
	plot(ladderize(tree), cex=0.35, label.offset=0.0005, tip.color=tipColors, y.lim=c(5,200))#, x.lim=c(-7,40))#, y.lim=c(130,3950))
	node.support(tree$posterior, mode="dots", col = "red", cex=0.3)
	node.support(tree$posterior, mode="dots", cutoff=0.95, col = "black", cex=0.55)
	node.support(tree$posterior, mode="numbers", digits=2,pos="below", col = "grey50", font=2, cex=0.3)

	add.scale.bar(x=par()$usr[1],y=1,lwd=1, cex=0.7)
	text(x=VAR,y=1,labels="expected substitutions per site", adj=0, cex=0.7)

	mtext(side=3, text=treeType, font=2, line=0, adj=0, cex=0.7)
	mtext(side=3, text="Plot #2 - Node support: black if PP >= 0.95; red if PP < 0.95 - Tips: red if ancient DNA", font=2, line= -1, adj=0, cex=0.7)
	#mtext(side=3, text="196 total tips (176 Microdipodops), 20 Dipodomys, excluding outgroup", font=3, line= -2, adj=0, cex=0.7)
	mtext(side=3, text="All Microdipodops, 20 Dipodomys, excluding outgroup", font=3, line= -2, adj=0, cex=0.7)

	#dev.off()


	## (3) MINUS DIPODOMYS 
	###
	#pdf(file=paste0("FIGURES/PHYLO/MCCtree_",cal,"_plot3_2022_noOutNoDipos_v2.pdf"), width=8.5, height=11, onefile=TRUE)

	tree<-drop.tip(treeToPlot, c("Perognathus_flavus", treeToPlot$tip.label[grep("Dipodomys",(treeToPlot$tip.label))]) )
	#correctPP_fullTree<-cbind(199:394,treeToPlot$posterior)
	#correctPP_tree2<-cbind(198:393,tree2$posterior)
	correctPP_tree3<-cbind(178:353,treeToPlot$posterior[21:196])

	#get tip colors
		tipColors<-rep("black",length(tree$tip.label))
		tipColors[match(aDNA_tips, tree$tip.label)] <- "red"

	#quartz(width=8.5, height=22)
	plot(ladderize(tree), cex=0.35, label.offset=0.00005, tip.color=tipColors, y.lim=c(5,180)) #x.lim=c(-7,40))#
	#node.support(tree$posterior, mode="dots", col = "red", cex=0.5)
	#node.support(tree$posterior, mode="dots", cutoff=0.95, col = "black", cex=0.55)
	#node.support(tree$posterior, mode="numbers", digits=2,pos="below", col = "black", font=2, cex=0.7)
	node.support(correctPP_tree3[,2], mode="dots", col = "red", cex=0.3)
	node.support(correctPP_tree3[,2], mode="dots", cutoff=0.95, col = "black", cex=0.55)
	node.support(correctPP_tree3[,2], mode="numbers", digits=2,pos="below", col = "grey50", font=2, cex=0.3)

	add.scale.bar(x=par()$usr[1],y=1,lwd=1, cex=0.7)
	text(x=VAR-0.02,y=1,labels="expected substitutions per site", adj=0, cex=0.7)

	mtext(side=3, text=treeType, font=2, line=0, adj=0, cex=0.7)
	mtext(side=3, text="Plot #3 - Node support: black if PP >= 0.95; red if PP < 0.95 - Tips: red if ancient DNA", font=2, line= -1, adj=0, cex=0.7)
	#mtext(side=3, text="176 Microdipodops only, excluding Dipodomys and outgroup", font=3, line= -2, adj=0, cex=0.7)
	mtext(side=3, text="All Microdipodops only, excluding Dipodomys and outgroup", font=3, line= -2, adj=0, cex=0.7)

	#dev.off()

	dev.off()






# CHRONOGRAMS:
#================

# Get the node support values + 95% HPDs per tree
#=================================================
	treeUni<-read_annotated(filename="final_reRuns_MB_dated-n-phylograms/infile_figtree.con_Micros_1Uni_labelsX_addSpace.tre", format="nexus")

	posterior<-c(); height_95_HPD_MIN<-c(); height_95_HPD_MAX<-c(); height_mean<-c()
	for(j in 1:length(treeUni$node.comment)){
		posterior[j]<-as.numeric(strsplit(strsplit(treeUni$node.comment[j],",")[[1]][1],"=")[[1]][2])
		height_95_HPD_MIN[j]<-as.numeric(strsplit(strsplit(treeUni$node.comment[j],",")[[1]][9],"=")[[1]][2])
		height_95_HPD_MAX[j]<-as.numeric(strsplit(treeUni$node.comment[j],",")[[1]][10])
		height_mean[j]<-as.numeric(strsplit(strsplit(treeUni$node.comment[j],",")[[1]][7],"=")[[1]][2])
	}
	treeUni$posterior <- posterior[(length(treeUni$tip.label)+1):length(treeUni$node.comment)]
	treeUni$'height_95%_HPD_MIN' <- height_95_HPD_MIN[(length(treeUni$tip.label)+1):length(treeUni$node.comment)]
	treeUni$'height_95%_HPD_MAX' <- height_95_HPD_MAX[(length(treeUni$tip.label)+1):length(treeUni$node.comment)]
	treeUni$height_mean <- height_mean[(length(treeUni$tip.label)+1):length(treeUni$node.comment)]

	treeUni_table<-cbind.data.frame(nodeNum=(length(treeUni$tip.label)+1):length(treeUni$node.comment),PP=treeUni$posterior, HPD95_min=treeUni$'height_95%_HPD_MIN', HPD95_max=treeUni$'height_95%_HPD_MAX', HPD95_mean=treeUni$height_mean)
	write.csv(treeUni_table, file="final_reRuns_MB_dated-n-phylograms/infile_figtree.con_Micros_1Uni_NodeSummaryTable.csv")
	pdf(file="final_reRuns_MB_dated-n-phylograms/infile_figtree.con_Micros_1Uni_NodeNumbersTree.pdf", width=40, height=20)
	plot(ladderize(treeUni), cex=0.3); nodelabels(cex=0.4)
	dev.off()

	treeExp<-read_annotated(filename="final_reRuns_MB_dated-n-phylograms/infile_figtree.con_Micros_1Exp_labelsX_addSpace.tre", format="nexus")
	posterior<-c(); height_95_HPD_MIN<-c(); height_95_HPD_MAX<-c(); height_mean<-c()
	for(j in 1:length(treeExp$node.comment)){
		posterior[j]<-as.numeric(strsplit(strsplit(treeExp$node.comment[j],",")[[1]][1],"=")[[1]][2])
		height_95_HPD_MIN[j]<-as.numeric(strsplit(strsplit(treeExp$node.comment[j],",")[[1]][9],"=")[[1]][2])
		height_95_HPD_MAX[j]<-as.numeric(strsplit(treeExp$node.comment[j],",")[[1]][10])
		height_mean[j]<-as.numeric(strsplit(strsplit(treeExp$node.comment[j],",")[[1]][7],"=")[[1]][2])
	}
	treeExp$posterior <- posterior[(length(treeExp$tip.label)+1):length(treeExp$node.comment)]
	treeExp$'height_95%_HPD_MIN' <- height_95_HPD_MIN[(length(treeExp$tip.label)+1):length(treeExp$node.comment)]
	treeExp$'height_95%_HPD_MAX' <- height_95_HPD_MAX[(length(treeExp$tip.label)+1):length(treeExp$node.comment)]
	treeExp$height_mean <- height_mean[(length(treeExp$tip.label)+1):length(treeExp$node.comment)]

	treeExp_table<-cbind.data.frame(nodeNum=(length(treeExp$tip.label)+1):length(treeExp$node.comment),PP=treeExp$posterior, HPD95_min=treeExp$'height_95%_HPD_MIN', HPD95_max=treeExp$'height_95%_HPD_MAX', HPD95_mean=treeExp$height_mean)
	write.csv(treeExp_table, file="final_reRuns_MB_dated-n-phylograms/infile_figtree.con_Micros_1Exp_NodeSummaryTable.csv")
	pdf(file="final_reRuns_MB_dated-n-phylograms/infile_figtree.con_Micros_1Exp_NodeNumbersTree_newTipNames.pdf", width=40, height=20)
	plot(ladderize(treeExp), cex=0.3); nodelabels(cex=0.4)
	dev.off()


# PLOT -- multi-page PDF with each of the iterations
#=====================	
	treeType<-"1 Exp calibration (soft maximum): birth-death, relaxed clock, MCC of 10k trees, mean node ages"
	calDetail<-"Dipodomyinae: 15.9 Ma min age, 20.6 Ma soft max age = 17.5 Ma mean of exponential prior"
	cal<-"dated-1Exp"
	treeToPlot<- treeExp

#	treeType<-"1 Uni calibration (hard maximum): birth-death, relaxed clock, MCC of 10k trees, mean node ages"
#	calDetail<-"Dipodomyinae: 15.9 Ma min age, 20.6 Ma hard max age"
#	cal<-"1Uni"
#	treeToPlot<-treeUni

	pdf(file=paste0("FIGURES/PHYLO/MCCtree_",cal,"_plot1-4_2022_ALL_v2.pdf"), width=8.5, height=11, onefile=TRUE)

	## (1) ALL SPECIES
	###
	#pdf(file=paste0("FIGURES/PHYLO/MCCtree_plot_Microdipodops_AgePPdots_",cal,"_ColxTips_2022_ALL.pdf"), width=20, height=40, onefile=TRUE)
	#pdf(file=paste0("FIGURES/PHYLO/MCCtree_",cal,"_plot1_2022_allTips_noPP_v2.pdf"), width=8.5, height=11, onefile=TRUE)

	tree<-treeToPlot
	#get tip colors
		sepNames<-do.call(rbind,strsplit(tree$tip.label,"_"))
		tipColors<-rep("black",length(tree$tip.label))
		tipColors[which(sepNames[,5]=="X")] <- "red"

	#quartz(width=8.5, height=22)
	plot(ladderize(tree), cex=0.3, label.offset=0.3, tip.color=tipColors, y.lim=c(5,195), x.lim=c(-3,39))

	HPDbars(tree, label="height_95%_HPD", broken=T, lwd=2, col=hsv(0.65,1,1,alpha=0.7))
	node.support(tree$posterior, mode="dots", col = "red", cex=0.3)
	node.support(tree$posterior, mode="dots", cutoff=0.95, col = "black", cex=0.55)
	node.support(branching.times(tree), mode="numbers", digits=1,pos="above", col = "black", font=2, cex=0.3)
	#node.support(tree$posterior, mode="numbers", digits=2,pos="below", col = "grey50", font=2, cex=0.5)

	data(gradstein04)
	axisGeo(GTS = gradstein04, unit = c("epoch","period"), cex = 0.5, ages=FALSE, gridty=3, gridcol="grey50")
	axisPhylo(cex.axis=0.5, pos=-5, mgp=c(0,0.2,0))

	mtext(side=3, text=treeType, font=2, line=1, adj=0, cex=0.7)
	#mtext(side=3, text=calDetail, font=1, line= 0, adj=0, cex=0.7)
	mtext(side=3, text="Plot #1 - Node support: black if PP >= 0.95; red if PP < 0.95 - Tips: red if ancient DNA", font=2, line= 0, adj=0, cex=0.7)
	mtext(side=3, text="197 total tips (176 Microdipodops), 20 Dipodomys, 1 Perognathus outgroup", font=3, line= -1, adj=0, cex=0.7)

	#dev.off()

	## (2) MINUS 1 OUTGROUP 
	###
	#pdf(file=paste0("FIGURES/PHYLO/MCCtree_",cal,"_plot2_2022_noOut_noPP_v2.pdf"), width=8.5, height=11, onefile=TRUE)

	#tree<-treeExp
	tree<-drop.tip2(treeToPlot,"Perognathus_flavus")
	#get tip colors
		sepNames<-do.call(rbind,strsplit(tree$tip.label,"_"))
		tipColors<-rep("black",length(tree$tip.label))
		tipColors[which(sepNames[,5]=="X")] <- "red"

	#quartz(width=8.5, height=22)
	plot(ladderize(tree), cex=0.35, label.offset=0.1, tip.color=tipColors, y.lim=c(5,195), x.lim=c(0,19))

	HPDbars(tree, label="height_95%_HPD", broken=T, lwd=2, col=hsv(0.65,1,1,alpha=0.7))
	node.support(tree$posterior, mode="dots", col = "red", cex=0.3)
	node.support(tree$posterior, mode="dots", cutoff=0.95, col = "black", cex=0.55)
	node.support(branching.times(tree), mode="numbers", digits=1,pos="above", col = "black", font=2, cex=0.3)
	#node.support(tree$posterior, mode="numbers", digits=2,pos="below", col = "grey50", font=2, cex=0.5)

	data(gradstein04)
	axisGeo(GTS = gradstein04, unit = c("epoch","period"), cex = 0.5, ages=FALSE, gridty=3, gridcol="grey50")
	axisPhylo(cex.axis=0.5, pos=-5, mgp=c(0,0.2,0))

	mtext(side=3, text=treeType, font=2, line=1, adj=0, cex=0.7)
	#mtext(side=3, text=calDetail, font=1, line= -1.5, adj=0, cex=0.7)
	mtext(side=3, text="Plot #2 - Node support: black if PP >= 0.95; red if PP < 0.95 - Tips: red if ancient DNA", font=2, line= 0, adj=0, cex=0.7)
	mtext(side=3, text="196 total tips (176 Microdipodops), 20 Dipodomys, excluding outgroup", font=3, line= -1, adj=0, cex=0.7)

	#dev.off()


	## (3) MINUS DIPODOMYS 
	###
	#pdf(file=paste0("FIGURES/PHYLO/MCCtree_",cal,"_plot3_2022_noOutNoDipos_noPP_v2.pdf"), width=8.5, height=11, onefile=TRUE)

	#tree<-treeExp
	tree<-drop.tip2(treeToPlot, c("Perognathus_flavus", treeExp$tip.label[grep("Dipodomys",(treeExp$tip.label))]) )
	#get tip colors
		sepNames<-do.call(rbind,strsplit(tree$tip.label,"_"))
		tipColors<-rep("black",length(tree$tip.label))
		tipColors[which(sepNames[,5]=="X")] <- "red"
	#quartz(width=8.5, height=22)
	plot(ladderize(tree), cex=0.37, label.offset=0.025, tip.color=tipColors, y.lim=c(5,170), x.lim=c(-0.4,5.5))

	HPDbars(tree, label="height_95%_HPD", broken=T, lwd=2, col=hsv(0.65,1,1,alpha=0.7))
	node.support(tree$posterior, mode="dots", col = "red", cex=0.3)
	node.support(tree$posterior, mode="dots", cutoff=0.95, col = "black", cex=0.55)
	node.support(branching.times(tree), mode="numbers", digits=1,pos="above", col = "black", font=2, cex=0.3)
	#node.support(tree$posterior, mode="numbers", digits=2,pos="below", col = "grey50", font=2, cex=0.5)

	data(gradstein04)
	axisGeo(GTS = gradstein04, unit = c("epoch","period"), cex = 0.5, ages=FALSE, gridty=3, gridcol="grey50")
	axisPhylo(cex.axis=0.5, pos=-5, mgp=c(0,0.2,0))

	mtext(side=3, text=treeType, font=2, line=1, adj=0, cex=0.7)
	#mtext(side=3, text=calDetail, font=1, line= -1.5, adj=0, cex=0.7)
	mtext(side=3, text="Plot #3 - Node support: black if PP >= 0.95; red if PP < 0.95 - Tips: red if ancient DNA", font=2, line= 0, adj=0, cex=0.7)
	mtext(side=3, text="176 Microdipodops only, excluding Dipodomys and outgroup", font=3, line= -1, adj=0, cex=0.7)

	#dev.off()


	## (4) Major clades + DIPODOMYS 
	###
	# keep only 1 rep per clade
	toKeep<-c(treeExp$tip.label[grep("Dipodomys",(treeExp$tip.label))], "Mina__MLZ_1782","Alamo__MSB_35536", "Owyhee__MLZ_2180", "FortRock__MLZ_2174","Minersville__MLZ_2071","WEureka__MLZ_2031")
	toDrop<-setdiff(treeExp$tip.label,toKeep)
	treeExp_simp <- drop.tip2(treeExp,toDrop)
		# relabel:
		oldLabel<-c("Mina__MLZ_1782","Alamo__MSB_35536", "Owyhee__MLZ_2180", "FortRock__MLZ_2174","Minersville__MLZ_2071","WEureka__MLZ_2031")
		newLabel<-c("M. (P.) pallidus", "M. (P.) ruficollaris", "M. (M.) megacephalus", "M. (M.) oregonus", "M. (M.) albiventer", "M. (M.) polionotus")
		for(i in 1: length(oldLabel)){
			tipNum<-match(oldLabel[i], treeExp_simp$tip.label)
			treeExp_simp$tip.label[tipNum]<-newLabel[i]
		}

	# PLOT simplified version...
	#pdf(file=paste0("FIGURES/PHYLO/MCCtree_",cal,"_plot4_2022_Micros-vs-Dipos_noPP_v2.pdf"), width=8.5, height=11, onefile=TRUE)

	tree<-treeExp_simp

	x1= -7 ; x2= 30
	y1= -15 ; y2= 30

	# plot dummy tree
	plot(ladderize(tree), cex=0.8, label.offset=0.4, tip.color="white", edge.width=2, y.lim=c(y1,y2), x.lim=c(x1,x2))
		# put rectangle
		root<-max(branching.times(tree))
		rect(xleft=root-2.6, ybottom=0, xright=root-0.8, ytop=27, col=grey(0.5, alpha=0.3), border=NA)

	par(new=TRUE)
	#quartz(width=8.5, height=22)
	plot(ladderize(tree), cex=0.8, label.offset=0.4, tip.color="black", edge.width=2, y.lim=c(y1,y2), x.lim=c(x1,x2))

	HPDbars(tree, label="height_95%_HPD", broken=T, lwd=3, col=hsv(0.65,1,1,alpha=0.7))
	node.support(tree$posterior, mode="dots", col = "red", cex=0.4)
	node.support(tree$posterior, mode="dots", cutoff=0.95, col = "black", cex=0.75)
	node.support(branching.times(tree), mode="numbers", digits=1,pos="above", col = "black", font=2, cex=0.8)
	#node.support(tree$posterior, mode="numbers", digits=2,pos="below", col = "grey50", font=2, cex=0.8)

	data(gradstein04)
	axisGeo(GTS = gradstein04, unit = c("epoch","period"), cex = 0.8, ages=FALSE, gridty=3, gridcol="grey50")
	axisPhylo(cex.axis=0.8, pos=-2)#, mgp=c(0,0.1,0))

	mtext(side=3, text=treeType, font=2, line=1, adj=0, cex=0.7)
	#mtext(side=3, text=calDetail, font=1, line= -1.5, adj=0, cex=0.7)
	mtext(side=3, text="Plot #4 - Node support: black if PP >= 0.95; red if PP < 0.95", font=2, line= 0, adj=0, cex=0.7)
	mtext(side=3, text="6 lineages of Microdipodops vs. 20 species of Dipodomys (gray box is 95% HPD of Microdipodops sister species)", font=3, line= -1, adj=0, cex=0.7)

	dev.off()







