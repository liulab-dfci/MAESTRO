## Marker gene list ##
srcdir = "../"
human_immune_simple <- list(B = c("MS4A1", "SLAMF7", "CD79A", "BLNK", "FCRL5"),
							  MacrophageM0 = c("PLA2G7", "PPBP", "QPCT", "SLAMF8", "SLC12A8", "TNFSF14"),
							  MacrophageM1 = c("PLA1A", "PTGIR", "RASSF4", "RSAD2", "SLAMF1", "SLC2A6","SOCS1","TNFAIP6","TNIP3","TRPM4"),
							  MacrophageM2 = c("CCL14","CCL23","CD4","CRYBB1","FRMD4A","GSTT1","HRH1","NPL","RENBP","WNT5B"),
							  CD14_Monocytes = c("LYZ", "CD14"),
							  FCGR3A_Monocytes = c("FCGR3A", "MS4A7"),
							  CD4T = c("IL7R", "TCF7", "SELL", "CCR7","GZMA","CCL5"),
							  Treg = c("FOXP3", "IL2RA", "TNFRSF9", "TIGIT", "CTLA4", "CCR8"),
							  CD8T = c("GZMA", "GZMB", "PRF1", "LEF1", "GZMH", "GZMK", "CCR7"),
							  Texh = c("CXCL13", "CTLA4", "TIGIT", "PDCD1", "HAVCR2", "ENTPD1",
							          "TNFRSF9", "CD27", "LAYN", "PHLDA1", "SNAP47","LAG3"),
							  NK = c("NKG7", "GNLY"),
							  Mast = c("CMA1", "MS4A2", "TPSAB1", "TPSB2"),
							  DC = c("CD40", "FCER1A", "CST3", "CD80", "CD83", "CCR7"),
					          Megakaryocytes = "PPBP",
					          Fibroblasts = c("FAP", "PDPN", "MMP2", "PDGFRA", "THY1", "MMP11","PDGFRL", "TGFB3", "COL1A2", "DCN", "COL3A1", "COL6A1"),
					          Myofibroblasts = c("ACTA2", "MCAM", "MYLK", "MYL9", "IL6", "PDGFA"),
					          Endothelial = c("PECAM1", "VWF", "ENG"))

human_stormal_simple <- list(Megakaryocytes = "PPBP",
					 		 Fibroblasts = c("FAP", "PDPN", "MMP2", "PDGFRA", "THY1", "MMP11","PDGFRL", "TGFB3", "COL1A2", "DCN", "COL3A1", "COL6A1"),
							 Myofibroblasts = c("ACTA2", "MCAM", "MYLK", "MYL9", "IL6", "PDGFA"),
					 		 Endothelial = c("PECAM1", "VWF", "ENG"))

TCIA <- read.delim(paste0(srcdir,"annotations/human_immune_TCIA.txt"))
TCIA_celltype = unique(TCIA$Celltype)
human_immune_TCIA = sapply(1:length(TCIA_celltype),function(i){return(as.character(TCIA[TCIA$Celltype == TCIA_celltype[i],2]))})
names(human_immune_TCIA) = TCIA_celltype
human_immune_TCIA <- c(human_immune_TCIA, human_immune_simple)

CIBERSORT <- read.delim(paste0(srcdir,"annotations/human_immune_CIBERSORT.txt"))
CIBERSORT_celltype = unique(CIBERSORT$CellType)
human_immune_CIBERSORT = sapply(1:length(CIBERSORT_celltype),function(i){return(as.character(CIBERSORT[CIBERSORT$CellType == CIBERSORT_celltype[i],2]))})
names(human_immune_CIBERSORT) = CIBERSORT_celltype
human_immune_CIBERSORT <- c(human_immune_CIBERSORT, human_immune_simple)

BioLegend = read.table(paste0(srcdir,"annotations/human_BioLegend.txt"),header = TRUE, sep = '\t',check.names = FALSE)
BioLegend_celltype = unique(BioLegend$Celltype)
human_immune_BioLegend = sapply(1:length(BioLegend_celltype),function(i){return(as.character(BioLegend[BioLegend$Celltype == BioLegend_celltype[i],2]))})
names(human_immune_BioLegend) = BioLegend_celltype
human_immune_BioLegend <- c(human_immune_BioLegend, human_immune_simple)

MouseBrain <- read.table(paste0(srcdir,"annotations/mouse_brain.txt"),header = TRUE,fill = TRUE, stringsAsFactors = FALSE,sep = "\t")
mouse_brain <- list(Ependymal = toupper(MouseBrain[MouseBrain[,1] != "",1]),
							Oligodendrocyte = toupper(MouseBrain[MouseBrain[,2] != "",2]),
							Microglia = toupper(MouseBrain[MouseBrain[,3] != "",3]),
							CA1Pyramidal = toupper(MouseBrain[MouseBrain[,4] != "",4]),
							Interneuron = toupper(MouseBrain[MouseBrain[,5] != "",5]),
							Endothelial = toupper(MouseBrain[MouseBrain[,6] != "",6]),
							S1Pyramidal = toupper(MouseBrain[MouseBrain[,7] != "",7]),
							Astrocyte = toupper(MouseBrain[MouseBrain[,8] != "",8]),
							Mural = toupper(MouseBrain[MouseBrain[,9] != "",9]))

TabulaMuris_droplet_table = read.table(paste0(srcdir,"annotations/TabulaMuris_droplet_markers_all.txt"),header = TRUE, sep = '\t',check.names = FALSE)
TabulaMuris_droplet_celltype = unique(TabulaMuris_droplet_table$cluster)
TabulaMuris_droplet = sapply(1:length(TabulaMuris_droplet_celltype),function(i){
  							  return(TabulaMuris_droplet_table[TabulaMuris_droplet_table$cluster == TabulaMuris_droplet_celltype[i],2])
  							  })
names(TabulaMuris_droplet) = TabulaMuris_droplet_celltype

TabulaMuris_facs_table = read.table(paste0(srcdir,"annotations/TabulaMuris_facs_markers_all.txt"),header = TRUE, sep = '\t',check.names = FALSE)
TabulaMuris_facs_celltype = unique(TabulaMuris_facs_table$cluster)
TabulaMuris_facs = sapply(1:length(TabulaMuris_facs_celltype),function(i){
  						   return(TabulaMuris_facs_table[TabulaMuris_facs_table$cluster == TabulaMuris_facs_celltype[i],2])
						   })
names(TabulaMuris_facs) = TabulaMuris_facs_celltype
				