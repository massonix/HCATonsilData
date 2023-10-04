# This script keeps track of all the annotation changes and versions,
# and saves a dictionary to be used in the updateAnnotations function

annotations_dictionary <- list(
  dict_20220215_to_20220619 = c(
    # MBC
    "MBC-like_nonproli" = "Precursor MBCs",
    "MBC-like_proli1" = "Precursor MBCs",
    "MBC-like_proli2"= "Reactivated proliferative MBCs",
    "MBC-like_proli3" = "Reactivated proliferative MBCs",
    "MBC-like_FCRL4+"= "Reactivated proliferative MBCs",
    "PC-precursors" = "PC committed Light Zone GCBC",
    "GC-Tfh-0X40" = "GC-Tfh-OX40",
    "non-GC-Tf-regs" = "Eff-Tregs-IL32",
    "GC-Tf-regs" = "Tfr",
    "IFN CD8 T" = "IFN+ CD8 T",
    
    # GCBC
    "LZ FDC" = "FDC",
    "DZ FDC" = "COL27A1+ FDC",
    "DZ_Sphase" = "DZ early Sphase",
    "DZ_Sphase_HistoneHigh" = "DZ late Sphase",
    "DZ_G2M_HistoneHigh" = "DZ early G2Mphase",
    "DZ_G2M_CCNBHigh"= "DZ late G2Mphase",
    "DZ-cell cycle exit" = "DZ cell cycle exit",
    "DZ-nonproliferative" = "DZ cell cycle exit",
    "DZ-nonproliferative_FOXP1hi"= "DZ non proliferative",
    "DZ/LZ" = "DZ_LZ transition",
    "DZ/LZ" = "DZ_LZ transition",
    "LZ" = "LZ",
    "LZ-BCL2A1 neg"= "LZ",
    "LZ-DZ-re-entry early commitment" = "LZ_DZ reentry commitment",
    "LZ-proliferative_BCL2A1pos" = "LZ proliferative",
    "LZ-proliferative_BCL2A1neg" = "LZ_DZ transition",
    
    # PC
    "class switch MBC" = "csMBC",
    "Neutrophil Granulocytes" = "Neutrophils",
    
    # CD8 T
    "CXCR6+ RM CD8 T" = "RM CD8 activated T",
    
    # myeloid
    "IL7R MMP12 macrophages" = "MMP Slancytes",
    "C1Q HLA macrophages" = "C1Q Slancytes",
    "SELENOP FUCA1 PTGDS macrophages" = "SELENOP Slancytes",
    "ITGAX ZEB2 macrophages" = "ITGAX Slancytes",
    "Mast cells" = "Mast"
  ),
  dict_20220619_to_20230508 = c(
    # Unconventional CD8 T cells
    "MAIT" = "MAIT/Vδ2+ γδ T",
    "TCRVδ+ gd T" = "non-Vδ2+ γδ T",
    "CD56+ gd T" = "ZNF683+ CD8 T",
    "Nksig CD8 T" = "EM CD8 T",
    
    # Slan-like
    "C1Q Slancytes" = "C1Q Slan-like",
    "ITGAX Slancytes" = "ITGAX Slan-like",
    "MMP Slancytes" = "MMP Slan-like",
    "SELENOP Slancytes" = "SELENOP Slan-like",
    
    # NK
    "CD16-CD56- NK" = "CD16-CD56dim NK"
  )
)

save(annotations_dictionary, file = here::here("data/annotations_dictionary.RData"))
