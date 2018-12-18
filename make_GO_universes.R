setwd("/home/pollardlab/Desktop/GO_annotations")
# This file imports gene ontology annotations from various sources, organizes them, and outputs .rds files with the GO universes to be used by GOStats in another R script.
# The GO annotations come from various sources for Tetrahymena, or are derived using BLASTp.


require(GO.db)
require(xlsx)
require(readxl)
require(gridExtra)


get_num_genes_added <- function(universe_to_compare, gene_ontology_org_df){
  # for making the table comparing the various sources of GO annotations to the geneontology.org source
  TTHERMs_from_test_df <- unique(as.character(universe_to_compare$frame.gene_id))
  TTHERMs_from_geneontology <- unique(as.character(gene_ontology_org_df$frame.gene_id))
  
  # get the TTHERMs unique to the first df
  TTHERMs_from_test_df_new <- TTHERMs_from_test_df[! TTHERMs_from_test_df %in% TTHERMs_from_geneontology]  # get the ones that are not in the gene_ontology_org_df
  
  result <- length(TTHERMs_from_test_df_new)
  return(result)
}


get_num_annotations_added <- function(universe_to_compare, gene_ontology_org_df){
  # for making the table comparing the various sources of GO annotations to the geneontology.org source
  annotations_from_test_df <- unique(paste(as.character(universe_to_compare$frame.go_id), as.character(universe_to_compare$frame.gene_id)))
  annotations_from_geneontology <- unique(paste(as.character(gene_ontology_org_df$frame.go_id), as.character(gene_ontology_org_df$frame.gene_id)))
  
  # get the TTHERMs unique to the first df
  annotations_from_test_df_new <- annotations_from_test_df[! annotations_from_test_df %in% annotations_from_geneontology]  # get the ones that are not in the gene_ontology_org_df
  
  result <- length(annotations_from_test_df_new)
  return(result)
}


make_BLAST_best_HSP_data.frame_from_fmt6_full_info <- function(file = "BLAST_output.txt", get_single_best_HSP = FALSE){
  # loads the blast output file and outputs a df with the hit HSP with the highest percent identity, this is a nuanced point as all HSPs have the same hit to GO term associations
  BLAST_csv <- read.csv(file, header = FALSE, sep = "\t")
  blast_df_total <- data.frame(Tet_gene = BLAST_csv$V1, frame.gene_id = BLAST_csv$V2,	e_val = BLAST_csv$V9, percent_identity = BLAST_csv$V6, alignment_length = BLAST_csv$V5, slen = BLAST_csv$V4, qlen = BLAST_csv$V3, nmatches = BLAST_csv$V10, ngaps = BLAST_csv$V8)
  
  # make new column names
  column_names_to_use <- c("Tet_gene", "frame.gene_id", "e_val", "percent_identity", "alignment_length", "slen", "qlen", "num_matches", "number_of_gaps")
  colnames(blast_df_total) <- column_names_to_use
  
  # get a single HSP for a given subject-querry match
  if (get_single_best_HSP == TRUE){
    blast_df_best_HSP <- data.frame(matrix(nrow = 0, ncol = length(column_names_to_use)))
    colnames(blast_df_best_HSP) <- column_names_to_use
    
    for (TTHERM in unique(blast_df_total$Tet_gene)){
      subsetted_blast_df <- blast_df_total[blast_df_total$Tet_gene == TTHERM, ]  # get the section of the blast table for that TTHERM <-> gene hit
      subsetted_blast_df <- subsetted_blast_df[with(subsetted_blast_df, order(-subsetted_blast_df$alignment_length)),]  # order the HSPs by alignment length
      subsetted_blast_df_best_HSP <- subsetted_blast_df[1,]
      blast_df_best_HSP <- rbind(blast_df_best_HSP, subsetted_blast_df_best_HSP)
    }
    return(blast_df_best_HSP)
  } else {
    return(blast_df_total)
  }
}


replace_common_name_with_TTHERM <- function(input_GO_universe_df){
  # some annotations have the common name instead of the TTHERM, so replace them with the TTHERM
  non_TTHERM_genes <- droplevels(unique(input_GO_universe_df$frame.gene_id[! grepl("TTHERM", input_GO_universe_df$frame.gene_id)]))
  for (non_TTHERM in non_TTHERM_genes){
    if (non_TTHERM %in% tetrahymena_gene_info_df$gene_name){  # if the gene is in the list of gene descriptions (not all are)
      correct_TTHERM <- as.character(droplevels(tetrahymena_gene_info_df[which(tetrahymena_gene_info_df$gene_name == non_TTHERM), 1]))  # get the correct TTHERM
      input_GO_universe_df[which(input_GO_universe_df$frame.gene_id == non_TTHERM) , 3] <- correct_TTHERM  # replace genes with the correct TTHERM
    }
  }
  return(input_GO_universe_df)
}


common_precessing_steps <- function(input_universe_df, rename_columns = FALSE){
  # performs some of the common processing steps required for each source of GO annotations
  if (rename_columns == TRUE){
    colnames(input_universe_df) <- c("frame.go_id",	"frame.Evidance"	,"frame.gene_id")
  }
  input_universe_df <- replace_common_name_with_TTHERM(input_universe_df)  # some annotations have the common name instead of the TTHERM, so replace them with the TTHERM
  input_universe_df <- input_universe_df[grepl("TTHERM", input_universe_df$frame.gene_id), ]  # remove GO annotations with genes that are not TTHERMs (e.g. legacy annotations: common names that are no longer associated with a TTHERM)
  input_universe_df <- input_universe_df[!duplicated(input_universe_df[c("frame.go_id","frame.gene_id")]), ] # remove go annotations that are redundant, not based on evidence
  return(input_universe_df)
}


##  LOAD ANNOTATIONS TO EXCLUDE  ##------------------------------------------------------------------------------------------------------------------------------------------------------------------
setwd("/home/pollardlab/Dropbox/Brian_Miller/handpicking_BLASTp-derived_GOterms_to_be_excluded")
GO_terms_to_exclude_df <- read_excel("BLASTp_derived_GO_terms_to_exclude_20181024.xlsx")
GO_terms_to_exclude <- as.character(GO_terms_to_exclude_df$GO_id)
setwd("/home/pollardlab/Desktop/GO_annotations")


##  TETRAHYMENA THERMOPHILA (GO annotations from various sources) ##-----------------------------------------------------------------------------------------------------------------------------------

# import tetrahymena gene description information
tet_gene_info_file <- "tetrahymena/tetrahymena_ttherm_descriptions_from_gff_and_gene_names_from_ciliate_org.txt"
tetrahymena_gene_info_df <-read.csv(tet_gene_info_file, header = TRUE, sep = ",")

# import GO annotations from http://geneontology.org (Amigo2)
geneOntology_org_filename <- "tetrahymena/geneontology_org_GO_annotations_downloaded20180920.txt"
geneonOntology_org_csv <- read.csv(geneOntology_org_filename, header=FALSE, sep="\t")
geneOntology_org_GO_universe <- data.frame(GOterms=geneonOntology_org_csv$V5, Evidence_Code=rep("TAS", nrow(geneonOntology_org_csv)), TTHERM = geneonOntology_org_csv$V3, stringsAsFactors = FALSE)
geneOntology_org_GO_universe <- common_precessing_steps(geneOntology_org_GO_universe, rename_columns = TRUE)
write.csv(geneOntology_org_GO_universe, file = "output_files/geneOntology_org_GO_universe.csv", row.names=FALSE)

# import GO annotations from Tetramine (http://tfgd.ihb.ac.cn/)
tetramine_filename <- "tetrahymena/Tetramine_FGD_GO_annotations_downloaded20180717.tsv"
tetramine_csv <- read.csv(tetramine_filename, header=TRUE, sep="\t")
tetramine_GO_universe <- data.frame(GOterms = tetramine_csv$Gene...GO.Annotation...Ontology.Term...Identifier, Evidence_Code = rep("NAS", nrow(tetramine_csv)), TTHERM=tetramine_csv$Gene...DB.identifier, stringsAsFactors = FALSE)  # no evidence code was supplied by Tetramine, thus all annotations were given the NAS evidence code
tetramine_GO_universe <- common_precessing_steps(tetramine_GO_universe, rename_columns = TRUE)
write.csv(tetramine_GO_universe, file = "output_files/tetramine_GO_universe.csv", row.names=FALSE)

# import GO annotations from http://ciliate.org
ciliate_org_filename <- "tet_GO_annotations_from_ciliate_org/ciliate_org_GO_universe.csv"
ciliate_org_GO_universe <- read.csv(ciliate_org_filename, header = TRUE, sep = ",")
ciliate_org_GO_universe <- common_precessing_steps(ciliate_org_GO_universe, rename_columns = TRUE)
ciliate_org_GO_universe[which(ciliate_org_GO_universe$frame.Evidance %in% c("(IC","IMp","iss","(ND")), 2]  = "TAS"  # replace the few truncated GO annotations
write.csv(ciliate_org_GO_universe, file = "output_files/ciliate_org_GO_universe.csv", row.names=FALSE)


##  ARABIDOPSIS THALIANA (protein sequences and GO annotations from https://www.arabidopsis.org/)  ##-----------------------------------------------------------------
# load BLAST output
arabidopsis_blast_df <- make_BLAST_best_HSP_data.frame_from_fmt6_full_info(file = "arabidopsis/TAIR10_Tet_best_hit_20180822.txt", get_single_best_HSP = TRUE)
arabidopsis_blast_df$frame.gene_id <- substr(as.character(arabidopsis_blast_df$frame.gene_id), 1,9)  # remove the splicoform information from the gene name

# plot BLAST stats
# hist(arabidopsis_blast_df$percent_identity, xlim = c(0,100))
# hist(arabidopsis_blast_df[arabidopsis_blast_df$alignment_length < 1500, ]$alignment_length, xlim = c(0,1500), breaks = seq(0,1500,50))  # plot the alignment length if <1500 (most were)
# hist(-log10(arabidopsis_blast_df$e_val))

# load annotations for organism that BLAST was run on
annotation_file <- "arabidopsis/ATH_GO_GOSLIM_downloaded20180813.txt"
annotations_csv <- read.csv(annotation_file, header = FALSE, sep = "\t")
annotations_csv <- annotations_csv[! is.na(annotations_csv$V7), ]  # drop some rows that did not have annotations
arabidopsis_annotations_df <- data.frame(frame.go_id = annotations_csv$V6, frame.Evidance = annotations_csv$V10,	frame.gene_id = annotations_csv$V1)

# make full BLAST-based GO tables
arabidopsis_annotations_df_with_BLAST_hits <- arabidopsis_annotations_df[arabidopsis_annotations_df$frame.gene_id %in% arabidopsis_blast_df$frame.gene_id, ]
arabidopsis_annotations_df_with_BLAST_hits_with_TTHERM <- merge.data.frame(arabidopsis_annotations_df_with_BLAST_hits, arabidopsis_blast_df[,c(1:5)], by = "frame.gene_id")  # merge blast results (only gene names) with the arabidopsis Go annotaitons
arabidopsis_BLAST_derived_GO_universe <- arabidopsis_annotations_df_with_BLAST_hits_with_TTHERM[ ,c(2:7)]  # remove the arabidopsis gene name column
names(arabidopsis_BLAST_derived_GO_universe)[names(arabidopsis_BLAST_derived_GO_universe) == "Tet_gene"] <- "frame.gene_id"  # rename TTHERM column to frame.gene_id
arabidopsis_BLAST_derived_GO_universe <- common_precessing_steps(arabidopsis_BLAST_derived_GO_universe)

# exclude the BLASTp-derived multicellular related GO terms, write to csv
arabidopsis_BLAST_derived_GO_universe <- arabidopsis_BLAST_derived_GO_universe[! arabidopsis_BLAST_derived_GO_universe$frame.go_id %in% GO_terms_to_exclude, ]
write.csv(arabidopsis_BLAST_derived_GO_universe[,c(1:3)], file = "output_files/arabidopsis_BLAST_derived_GO_universe.csv", row.names=FALSE)


##  DROSOPHILA MELANGASTER (protein sequences from http://flybase.org, GO annotations from http://geneontology.org)  ##-------------------------------------------------
# load BLAST output
drosophila_blast_df <- make_BLAST_best_HSP_data.frame_from_fmt6_full_info(file = "drosophila_melanogaster/dmel_Tet_best_hit_20180923.txt", get_single_best_HSP = TRUE)

# plot BLAST stats
# hist(drosophila_blast_df$percent_identity, xlim = c(0,100))
# hist(drosophila_blast_df$alignment_length, xlim = c(0,2000), breaks = seq(0,1500,50))
# hist(drosophila_blast_df[drosophila_blast_df$alignment_length < 1500, ]$alignment_length, xlim = c(0,1500), breaks = seq(0,1500,50))  # plot the alignment length if <1500 (most were)
# hist(-log10(drosophila_blast_df$e_val))

# load annotations for organism that BLAST was run on
annotation_file <- "drosophila_melanogaster/fb_geneontology_org_20180815.gaf"
drosophila_annotations_csv <- read.csv(annotation_file, sep = "\t", skip = 16, header = FALSE)
drosophila_annotations_df <- data.frame(frame.go_id = drosophila_annotations_csv$V5, frame.Evidance = drosophila_annotations_csv$V7,	frame.gene_id = drosophila_annotations_csv$V2)

# make full BLAST-based GO tables
drosophila_annotations_df_with_BLAST_hits <- drosophila_annotations_df[drosophila_annotations_df$frame.gene_id %in% drosophila_blast_df$frame.gene_id, ]
drosophila_annotations_df_with_BLAST_hits_with_TTHERM <- merge.data.frame(drosophila_annotations_df_with_BLAST_hits, drosophila_blast_df[,c(1:5)], by = "frame.gene_id")  # merge blast results (only gene names) with the drosophila Go annotaitons
drosophila_BLAST_derived_GO_universe <- drosophila_annotations_df_with_BLAST_hits_with_TTHERM[ ,c(2:7)]  # remove the drosophila gene name column
names(drosophila_BLAST_derived_GO_universe)[names(drosophila_BLAST_derived_GO_universe) == "Tet_gene"] <- "frame.gene_id"  # rename TTHERM column to frame.gene_id
drosophila_BLAST_derived_GO_universe <- common_precessing_steps(drosophila_BLAST_derived_GO_universe)

# exclude the BLASTp-derived multicellular related GO terms, write to csv
drosophila_BLAST_derived_GO_universe <- drosophila_BLAST_derived_GO_universe[!drosophila_BLAST_derived_GO_universe$frame.go_id %in% GO_terms_to_exclude, ]
write.csv(drosophila_BLAST_derived_GO_universe[,c(1:3)], file = "output_files/drosophila_BLAST_derived_GO_universe.csv", row.names=FALSE)


##  PLASMODIUM FALCIPARUM (protein sequences from https://www.uniprot.org/)  ##------------------------------------------------------------------------------------------------
# load BLAST output
plasmodium_blast_df <- make_BLAST_best_HSP_data.frame_from_fmt6_full_info(file = "plasmodium_falciparum/plasmodium_Tet_best_hit_20180822.txt")

# plot BLAST stats
# hist(plasmodium_blast_df$percent_identity, xlim = c(0,100))
# hist(-log10(plasmodium_blast_df$e_val))

# load annotations for organism that BLAST was run on
annotation_file <- "plasmodium_falciparum/genedb_pfalciparum_valid.gaf"
plasmodium_annotations_csv <- read.csv(annotation_file, sep = "\t", skip = 5, header = FALSE)
plasmodium_annotations_df <-data.frame(frame.go_id = plasmodium_annotations_csv$V5, frame.Evidance = plasmodium_annotations_csv$V7,	frame.gene_id = plasmodium_annotations_csv$V2)
  
# make full BLAST-based GO tables
plasmodium_annotations_df_with_BLAST_hits <- plasmodium_annotations_df[plasmodium_annotations_df$frame.gene_id %in% plasmodium_blast_df$frame.gene_id, ]
plasmodium_annotations_df_with_BLAST_hits_with_TTHERM <- merge.data.frame(plasmodium_annotations_df_with_BLAST_hits, plasmodium_blast_df[,c(1:5)], by = "frame.gene_id")  # merge blast results (only gene names) with the plasmodium Go annotaitons
plasmodium_BLAST_derived_GO_universe <- plasmodium_annotations_df_with_BLAST_hits_with_TTHERM[ ,c(2:7)]  # remove the plasmodium gene name column
names(plasmodium_BLAST_derived_GO_universe)[names(plasmodium_BLAST_derived_GO_universe) == "Tet_gene"] <- "frame.gene_id"  # rename TTHERM column to frame.gene_id
plasmodium_BLAST_derived_GO_universe <- common_precessing_steps(plasmodium_BLAST_derived_GO_universe)

# exclude the BLASTp-derived multicellular related GO terms, write to csv
plasmodium_BLAST_derived_GO_universe <- plasmodium_BLAST_derived_GO_universe[! plasmodium_BLAST_derived_GO_universe$frame.go_id %in% GO_terms_to_exclude, ]
write.csv(plasmodium_BLAST_derived_GO_universe[,c(1:3)], file = "output_files/plasmodium_BLAST_derived_GO_universe.csv", row.names=FALSE)

# make txt file with the genes to supply to uniprot to get the protein sequences
# plasmodium_gene_names_for_uniprot <- as.character(plasmodium_annotations_df$frame.gene_id)
# fileConn<-file("plasmodium_gene_names_for_uniprot.txt")
# writeLines(plasmodium_gene_names_for_uniprot, fileConn)
# close(fileConn)


##  MERGE THE GO ANNOTATION FILES FOR ALL ORGANISMS (including Tetrahymena) and make a combined csv file  ##----------------------------------------------------------------------------------------------------------
all_sources_combined_GO_universe_noGeneOntology <- rbind(ciliate_org_GO_universe, tetramine_GO_universe)
all_sources_combined_GO_universe_yesGeneOntology <- rbind(all_sources_combined_GO_universe_noGeneOntology, geneOntology_org_GO_universe)
all_sources_combined_GO_universe_yesGeneOntology <- all_sources_combined_GO_universe_yesGeneOntology[!duplicated(all_sources_combined_GO_universe_yesGeneOntology[c("frame.go_id","frame.gene_id")]), ]  # remove redundancies
all_sources_combined_GO_universe_yesGeneOntology <- all_sources_combined_GO_universe_yesGeneOntology[! all_sources_combined_GO_universe_yesGeneOntology$frame.go_id %in% GO_terms_to_exclude, ]

# merge full BLAST GO df with the existing Tet GO anntoations
all_sources_combined_GO_universe_noGeneOntology <- rbind(all_sources_combined_GO_universe_noGeneOntology, arabidopsis_BLAST_derived_GO_universe[,c(1:3)])
all_sources_combined_GO_universe_noGeneOntology <- rbind(all_sources_combined_GO_universe_noGeneOntology, drosophila_BLAST_derived_GO_universe[,c(1:3)])
all_sources_combined_GO_universe_noGeneOntology <- rbind(all_sources_combined_GO_universe_noGeneOntology, plasmodium_BLAST_derived_GO_universe[,c(1:3)])

# remove redundant annotations
all_sources_combined_GO_universe_noGeneOntology <- all_sources_combined_GO_universe_noGeneOntology[!duplicated(all_sources_combined_GO_universe_noGeneOntology[c("frame.go_id","frame.gene_id")]), ]  # remove redundancies

# add the geneontology_org universe and remove redundancies / multicellular organisms
all_sources_combined_GO_universe <- rbind(all_sources_combined_GO_universe_noGeneOntology, geneOntology_org_GO_universe)  # add the GeneOntology.org annotations
all_sources_combined_GO_universe <- all_sources_combined_GO_universe[!duplicated(all_sources_combined_GO_universe[c("frame.go_id","frame.gene_id")]), ]  # remove redundancies
all_sources_combined_GO_universe <- all_sources_combined_GO_universe[! all_sources_combined_GO_universe$frame.go_id %in% GO_terms_to_exclude, ]

# write to csv
write.csv(all_sources_combined_GO_universe, file = "output_files/all_sources_combined_GO_universe.csv", row.names=FALSE)
saveRDS(all_sources_combined_GO_universe, file = "output_files/all_sources_combined_GO_universe_withBLAST.rds")
saveRDS(all_sources_combined_GO_universe_yesGeneOntology, file = "output_files/all_Tet_sources_combined_GO_universe.rds")

##  MAKING QC PLOTS USING THE BLASTp-DERIVED ANNOTATIONS  ##------------------------------------------------------------------------------------------------------------------------------------
# all_sources_combined_GO_universe_for_QC <- rbind(arabidopsis_BLAST_derived_GO_universe, drosophila_BLAST_derived_GO_universe)
# all_sources_combined_GO_universe_for_QC <- rbind(all_sources_combined_GO_universe_for_QC, plasmodium_BLAST_derived_GO_universe)

# plot BLAST stats, these are the combined BLAST annotations, the redundancies have not been removed
# hist(all_sources_combined_GO_universe_for_QC$percent_identity, xlim = c(0,100), breaks = seq(0,100,1))  # percent identity
# hist(all_sources_combined_GO_universe_for_QC[all_sources_combined_GO_universe_for_QC$alignment_length < 2000, ]$alignment_length, xlim = c(0,2000), breaks = seq(0,2000,50))  # plot the alignment length if <2000 (most were)
# hist(-log10(all_sources_combined_GO_universe_for_QC$e_val), breaks = seq(0,200,5))  # eval


## MAKE COMPARISON TABLE OF GO ANNOTATIONS  ##--------------------------------------------------------------------------------------------------------------------------------------------------------
# table 2
output_table2 <- data.frame(matrix(nrow = 5, ncol = 3))
colnames(output_table2) = c("Source Added to geneontology.org", "Annotations Added to geneontology.org Annotations", "Genes Added to geneontology.org Annotations")
output_table2[,1] = c("ciliate.org", "Tetramine", "BLASTp using A. thaliana*", "BLASTp using D. melanogaster*", "BLASTp using P. falciparum*")

# populate the table
output_table2[1,3] <- get_num_genes_added(ciliate_org_GO_universe, geneOntology_org_GO_universe)
output_table2[2,3] <- get_num_genes_added(tetramine_GO_universe, geneOntology_org_GO_universe)
output_table2[3,3] <- get_num_genes_added(arabidopsis_BLAST_derived_GO_universe[,c(1:3)], geneOntology_org_GO_universe)
output_table2[4,3] <- get_num_genes_added(drosophila_BLAST_derived_GO_universe[,c(1:3)], geneOntology_org_GO_universe)
output_table2[5,3] <- get_num_genes_added(plasmodium_BLAST_derived_GO_universe[,c(1:3)], geneOntology_org_GO_universe)

output_table2[1,2] <- get_num_annotations_added(ciliate_org_GO_universe, geneOntology_org_GO_universe)
output_table2[2,2] <- get_num_annotations_added(tetramine_GO_universe, geneOntology_org_GO_universe)
output_table2[3,2] <- get_num_annotations_added(arabidopsis_BLAST_derived_GO_universe[,c(1:3)], geneOntology_org_GO_universe)
output_table2[4,2] <- get_num_annotations_added(drosophila_BLAST_derived_GO_universe[,c(1:3)], geneOntology_org_GO_universe)
output_table2[5,2] <- get_num_annotations_added(plasmodium_BLAST_derived_GO_universe[,c(1:3)], geneOntology_org_GO_universe)

pdf("output_table2.pdf", height=5, width=12)
grid.table(output_table2,  rows = NULL)
dev.off()

# table 1
output_table1 <- data.frame(matrix(nrow = 3, ncol = 3))
colnames(output_table1) = c("Source of Annotations", "Number of GO Annotations", "Number of Genes Represented")
output_table1[,1] = c("geneontology.org", "All Tetrahymena Sources","All Sources Combined*")
output_table1[1,2] <- nrow(geneOntology_org_GO_universe)
output_table1[2,2] <- nrow(all_sources_combined_GO_universe_yesGeneOntology)
output_table1[3,2] <- nrow(all_sources_combined_GO_universe) 
output_table1[1,3] <- length(unique(geneOntology_org_GO_universe$frame.gene_id))
output_table1[2,3] <- length(unique(all_sources_combined_GO_universe_yesGeneOntology$frame.gene_id))
output_table1[3,3] <- length(unique(all_sources_combined_GO_universe$frame.gene_id))
pdf("output_table1.pdf", height=5, width=12)
grid.table(output_table1,  rows = NULL)
dev.off()
