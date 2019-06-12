## Marker gene list ##
srcdir = "../src/"
TCIA <- read.delim(paste0(srcdir,"annotations/human_immune_TCIA.txt"))
markers.TCIA <- list(CD8Ta = as.character(TCIA[which(TCIA[,1]=="ActivatedCD8Tcell"),2]),
					 CD8Tem = as.character(TCIA[which(TCIA[,1]=="EffectormemeoryCD8Tcell"),2]),
					 CD8Tcm = as.character(TCIA[which(TCIA[,1]=="CentralmemeoryCD8Tcell"),2]),
					 NKT = as.character(TCIA[which(TCIA[,1]=="NaturalkillerTcell"),2]),
					 CD4Ta = as.character(TCIA[which(TCIA[,1]=="ActivatedCD4Tcell"),2]),
					 CD4Tem = as.character(TCIA[which(TCIA[,1]=="EffectormemeoryCD4Tcell"),2]),
					 CD4Tcm = as.character(TCIA[which(TCIA[,1]=="CentralmemeoryCD4Tcell"),2]),				
					 Th1 = as.character(TCIA[which(TCIA[,1]=="Type1Thelpercell"),2]),
					 Th2 = as.character(TCIA[which(TCIA[,1]=="Type2Thelpercell"),2]),
					 Th17 = as.character(TCIA[which(TCIA[,1]=="Type17Thelpercell"),2]),
					 Tfh = as.character(TCIA[which(TCIA[,1]=="Tfollicularhelpercell"),2]),
					 Treg = as.character(TCIA[which(TCIA[,1]=="RegulatoryTcell"),2]),
					 Tgd = as.character(TCIA[which(TCIA[,1]=="GammadeltaTcell"),2]),
					 Bactivate = as.character(TCIA[which(TCIA[,1]=="ActivatedBcell"),2]),
					 Bmemory = as.character(TCIA[which(TCIA[,1]=="MemoryBcell"),2]),
					 Bnaive = as.character(TCIA[which(TCIA[,1]=="ImmatureBcell"),2]),
					 DCactivate = as.character(TCIA[which(TCIA[,1]=="Activateddendriticcell"),2]),
					 DCnaive = as.character(TCIA[which(TCIA[,1]=="Immaturedendriticcell"),2]),
					 DCPlasmacytoid = as.character(TCIA[which(TCIA[,1]=="Plasmacytoiddendriticcell"),2]),
					 Eosinophil = as.character(TCIA[which(TCIA[,1]=="Eosinophil"),2]),
					 Neutrophil = as.character(TCIA[which(TCIA[,1]=="Neutrophil"),2]),
					 NK = as.character(TCIA[which(TCIA[,1]=="Naturalkillercell"),2]),
					 NKCD56dim = as.character(TCIA[which(TCIA[,1]=="CD56dimnaturalkillercell"),2]),
					 NKCD56bright = as.character(TCIA[which(TCIA[,1]=="CD56brightnaturalkillercell"),2]),
					 Monocyte = as.character(TCIA[which(TCIA[,1]=="Monocyte"),2]),
					 Macrophage = as.character(TCIA[which(TCIA[,1]=="Macrophage"),2]),
					 Mast = as.character(TCIA[which(TCIA[,1]=="Mast cell"),2]),
					 MDSC = as.character(TCIA[which(TCIA[,1]=="MDSC"),2]),
					 Megakaryocytes = "PPBP",
					 CAF = c("FAP", "PDPN", "MMP2", "PDGFRA", "THY1", "MMP11","PDGFRL", "TGFB3"),
					 Fibroblasts = c("FAP", "PDPN", "COL1A2", "DCN", "COL3A1", "COL6A1"),
					 Myofibroblasts = c("ACTA2", "MCAM", "MYLK", "MYL9", "IL6", "PDGFA"),
					 Myocyte = c("ACTA1", "ACTN2", "MYL2", "MYH2"),
					 Endothelial = c("PECAM1", "VWF", "ENG"))

CIBERSORT <- read.delim(paste0(srcdir,"annotations/human_immune_CIBERSORT.txt"))
markers.CIBERSORT <- list(Bnaive = as.character(CIBERSORT[which(CIBERSORT[,1]=="NaiveBcell"),2]),
					 	  Bmemory = as.character(CIBERSORT[which(CIBERSORT[,1]=="MemoryBcell"),2]),
					 	  Bactivate = as.character(CIBERSORT[which(CIBERSORT[,1]=="ActivatedBcell"),2]),
					 	  CD8T = as.character(CIBERSORT[which(CIBERSORT[,1]=="CD8Tcell"),2]),
					 	  CD4Tnaive = as.character(CIBERSORT[which(CIBERSORT[,1]=="NaiveCD4Tcell"),2]),
					 	  CD4Tmresting = as.character(CIBERSORT[which(CIBERSORT[,1]=="RestingmemoryCD4Tcell"),2]),
					 	  CD4Tmactivate = as.character(CIBERSORT[which(CIBERSORT[,1]=="ActivatememoryCD4Tcell"),2]),				
					 	  Tfh = as.character(CIBERSORT[which(CIBERSORT[,1]=="Tfollicularhelpercell"),2]),
					 	  Treg = as.character(CIBERSORT[which(CIBERSORT[,1]=="RegulatoryTcell"),2]),
					 	  NKresting = as.character(CIBERSORT[which(CIBERSORT[,1]=="Restingnaturalkillercell"),2]),
					 	  NKactivate = as.character(CIBERSORT[which(CIBERSORT[,1]=="Activatenaturalkillercell"),2]),
					 	  Monocyte = as.character(CIBERSORT[which(CIBERSORT[,1]=="Monocyte"),2]),
					 	  MacrophageM0 = as.character(CIBERSORT[which(CIBERSORT[,1]=="MacrophageM0"),2]),
					 	  MacrophageM1 = as.character(CIBERSORT[which(CIBERSORT[,1]=="MacrophageM1"),2]),
					 	  MacrophageM2 = as.character(CIBERSORT[which(CIBERSORT[,1]=="MacrophageM2"),2]),
					 	  DCresting = as.character(CIBERSORT[which(CIBERSORT[,1]=="Restingdendriticcell"),2]),
					 	  DCactivate = as.character(CIBERSORT[which(CIBERSORT[,1]=="Activatedendriticcell"),2]),
					 	  Mastresting = as.character(CIBERSORT[which(CIBERSORT[,1]=="RestingMastcells"),2]),
					 	  Mastactivate = as.character(CIBERSORT[which(CIBERSORT[,1]=="ActivateMastcells"),2]),
					 	  Eosinophil = as.character(CIBERSORT[which(CIBERSORT[,1]=="Eosinophil"),2]),
					 	  Neutrophil = as.character(CIBERSORT[which(CIBERSORT[,1]=="Neutrophil"),2]),
					 	  Megakaryocytes = "PPBP",
					 	  CAF = c("FAP", "PDPN", "MMP2", "PDGFRA", "THY1", "MMP11","PDGFRL", "TGFB3"),
					 	  Fibroblasts = c("FAP", "PDPN", "COL1A2", "DCN", "COL3A1", "COL6A1"),
					 	  Myofibroblasts = c("ACTA2", "MCAM", "MYLK", "MYL9", "IL6", "PDGFA"),
					 	  Myocyte = c("ACTA1", "ACTN2", "MYL2", "MYH2"),
					 	  Endothelial = c("PECAM1", "VWF", "ENG"),
					      Texh = c("CXCL13", "CTLA4", "TIGIT", "PDCD1", "HAVCR2", "ENTPD1",
					           "TNFRSF9", "CD27", "LAYN", "PHLDA1", "SNAP47","LAG3"),
					      MAIT = c("SLC4A10", "ZBTB16", "RORC"),
					      Treg = c("FOXP3", "IL2RA", "TNFRSF9", "TIGIT", "CTLA4", "CCR8")
)

markers.immune.simple <- list(B = c("MS4A1", "SLAMF7", "CD79A", "BLNK", "FCRL5"),
							  MacrophageM0 = c("PLA2G7", "PPBP", "QPCT", "SLAMF8", "SLC12A8", "TNFSF14"),
							  MacrophageM1 = c("PLA1A", "PTGIR", "RASSF4", "RSAD2", "SLAMF1", "SLC2A6","SOCS1","TNFAIP6","TNIP3","TRPM4"),
							  MacrophageM2 = c("CCL14","CCL23","CD4","CRYBB1","FRMD4A","GSTT1","HRH1","NPL","RENBP","WNT5B"),
							  CD14.Monocytes = c("LYZ", "CD14"),
							  FCGR3A.Monocytes = c("FCGR3A", "MS4A7"),
							  CD4T = c("IL7R", "TCF7", "SELL", "CCR7","GZMA","CCL5"),
							  Treg = c("FOXP3", "IL2RA", "TNFRSF9", "TIGIT", "CTLA4", "CCR8"),
							  CD8T = c("GZMA", "GZMB", "PRF1", "LEF1", "GZMH", "GZMK", "CCR7"),
							  Texh = c("CXCL13", "CTLA4", "TIGIT", "PDCD1", "HAVCR2", "ENTPD1",
							          "TNFRSF9", "CD27", "LAYN", "PHLDA1", "SNAP47","LAG3"),
							  NK = c("NKG7", "GNLY"),
							  Mast = c("CMA1", "MS4A2", "TPSAB1", "TPSB2"),
							  DC = c("CD40", "FCER1A", "CST3", "CD80", "CD83", "CCR7"),
							  MAIT = c("SLC4A10", "ZBTB16", "RORC"),
							  Megakaryocytes = "PPBP",
							  CAF = c("FAP", "PDPN", "MMP2", "PDGFRA", "THY1", "MMP11","PDGFRL", "TGFB3"),
							  Fibroblasts = c("FAP", "PDPN", "COL1A2", "DCN", "COL3A1", "COL6A1"),
							  Myofibroblasts = c("ACTA2", "MCAM", "MYLK", "MYL9", "IL6", "PDGFA"),
							  Myocyte = c("ACTA1", "ACTN2", "MYL2", "MYH2"),
							  Endothelial = c("PECAM1", "VWF", "ENG")
)

MouseBrain <- read.table(paste0(srcdir,"annotations/mouse_brain.txt"),header = TRUE,fill = TRUE, stringsAsFactors = FALSE,sep = "\t")
markers.brain.adult <- list(Ependymal = toupper(MouseBrain[MouseBrain[,1] != "",1]),
							Oligodendrocyte = toupper(MouseBrain[MouseBrain[,2] != "",2]),
							Microglia = toupper(MouseBrain[MouseBrain[,3] != "",3]),
							CA1Pyramidal = toupper(MouseBrain[MouseBrain[,4] != "",4]),
							Interneuron = toupper(MouseBrain[MouseBrain[,5] != "",5]),
							Endothelial = toupper(MouseBrain[MouseBrain[,6] != "",6]),
							S1Pyramidal = toupper(MouseBrain[MouseBrain[,7] != "",7]),
							Astrocyte = toupper(MouseBrain[MouseBrain[,8] != "",8]),
							Mural = toupper(MouseBrain[MouseBrain[,9] != "",9]))

TabulaMuris_droplet = read.table(paste0(srcdir,"annotations/TabulaMuris_droplet_markers_all.txt"),header = TRUE, sep = '\t',check.names = FALSE)
TabulaMuris_facs = read.table(paste0(srcdir,"annotations/TabulaMuris_facs_markers_all.txt"),header = TRUE, sep = '\t',check.names = FALSE)
TabulaMuris_droplet_celltype = unique(TabulaMuris_droplet$cluster)
TabulaMuris.droplet.markers = sapply(1:length(TabulaMuris_droplet_celltype),function(i){
  							  return(TabulaMuris_droplet[TabulaMuris_droplet$cluster == TabulaMuris_droplet_celltype[i],2])
  							  })
names(TabulaMuris.droplet.markers) = TabulaMuris_droplet_celltype

TabulaMuris_facs_celltype = unique(TabulaMuris_facs$cluster)
TabulaMuris.facs.markers = sapply(1:length(TabulaMuris_facs_celltype),function(i){
  						   return(TabulaMuris_facs[TabulaMuris_facs$cluster == TabulaMuris_facs_celltype[i],2])
						   })
names(TabulaMuris.facs.markers) = TabulaMuris_facs_celltype

BioLegend = read.table(paste0(srcdir,"annotations/human_immune_BioLegend.txt"),header = TRUE, sep = '\t',check.names = FALSE)
BioLegend_celltype = unique(BioLegend$CellType)
markers.BioLegend = sapply(1:length(BioLegend_celltype),function(i){
  return(BioLegend[BioLegend$CellType == BioLegend_celltype[i],2])
})
names(markers.BioLegend) = BioLegend_celltype

# markers.Gabriela <- list(B = as.character(Gabriela[which(Gabriela[,1]=="Gabriela"),2]),
# 					 	 CD8T = as.character(Gabriela[which(Gabriela[,1]=="CD8Tcell"),2]),
# 					 	 CytoT = as.character(Gabriela[which(Gabriela[,1]=="Cytotoxiscell"),2]),
# 					 	 DC = as.character(Gabriela[which(Gabriela[,1]=="DC"),2]),
# 					 	 Eosinophil = as.character(Gabriela[which(Gabriela[,1]=="Eosinophils"),2]),
# 					 	 Macrophage = as.character(Gabriela[which(Gabriela[,1]=="Macrophage"),2]),
# 					 	 Mast = as.character(Gabriela[which(Gabriela[,1]=="Mastcell"),2]),				
# 					 	 NKCD56bright = as.character(Gabriela[which(Gabriela[,1]=="CD56brightnaturalkillercell"),2]),
# 					 	 NKCD56dim = as.character(Gabriela[which(Gabriela[,1]=="CD56dimnaturalkillercell"),2]),
# 					 	 NK = as.character(Gabriela[which(Gabriela[,1]=="Naturalkillercell"),2]),
# 					 	 Neutrophil = as.character(Gabriela[which(Gabriela[,1]=="Neutrophils"),2]),
# 					 	 Tfh = as.character(Gabriela[which(Gabriela[,1]=="Tfollicularhelpercell"),2]),
# 					 	 T = as.character(Gabriela[which(Gabriela[,1]=="Tcell"),2]),
# 					 	 Tcm = as.character(Gabriela[which(Gabriela[,1]=="Tcm"),2]),
# 					 	 Tem = as.character(Gabriela[which(Gabriela[,1]=="Tem"),2]),
# 					 	 Tgd = as.character(Gabriela[which(Gabriela[,1]=="GammadeltaTcell"),2]),
# 					 	 Th17 = as.character(Gabriela[which(Gabriela[,1]=="Type17Thelpercell"),2]),
# 					 	 Th1 = as.character(Gabriela[which(Gabriela[,1]=="Type1Thelpercell"),2]),
# 					 	 Th2 = as.character(Gabriela[which(Gabriela[,1]=="Type2Thelpercell"),2]),
# 					 	 Th = as.character(Gabriela[which(Gabriela[,1]=="Thelpercell"),2]),
# 					 	 Treg = as.character(Gabriela[which(Gabriela[,1]=="RegulatoryTcell"),2]),
# 					 	 aDC = as.character(Gabriela[which(Gabriela[,1]=="aDC"),2]),
# 					 	 iDC = as.character(Gabriela[which(Gabriela[,1]=="iDC"),2]),
# 					 	 pDC = as.character(Gabriela[which(Gabriela[,1]=="Plasmacytoiddendriticcell"),2]),
# 					 	 MAIT = c("SLC4A10", "ZBTB16", "RORC"),
# 					 	 Megakaryocytes = "PPBP",
# 					 	 CAF = c("FAP", "PDPN", "MMP2", "PDGFRA", "THY1", "MMP11","PDGFRL", "TGFB3"),
# 					 	 Fibroblasts = c("FAP", "PDPN", "COL1A2", "DCN", "COL3A1", "COL6A1"),
# 					 	 Myofibroblasts = c("ACTA2", "MCAM", "MYLK", "MYL9", "IL6", "PDGFA"),
# 					 	 Myocyte = c("ACTA1", "ACTN2", "MYL2", "MYH2"),
# 					 	 Endothelial = c("PECAM1", "VWF", "ENG"),
# 					 	 Tcellexhausted = c("CXCL13", "CTLA4", "TIGIT", "PDCD1", "HAVCR2",
# 									    "ENTPD1", "TNFRSF9", "CD27", "LAYN", "PHLDA1", "SNAP47",
# 									    "LAG3", "PD1"))

				