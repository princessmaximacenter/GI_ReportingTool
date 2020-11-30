# R script containing functions for processing information of candidate 
# genetic interaction.
#
# Author(s): Denise Kersjes
# Date of creation:  21  September 2020
# Date of last edit: 30 November 2020


### Load required files
vep.file <- 'data/vep_out_muts_target_dkfz.txt'
muts.file <- 'data/muts_cand_target_dkfz.txt'
cand.file <- 'data/cand_target_dkfz.txt'
vep.df <- utils::read.delim(file = vep.file, header = T, 
                            stringsAsFactors = F, sep = "\t")
muts.df <- utils::read.delim(file = muts.file, header = T, 
                             stringsAsFactors = F, sep = "\t")
cand.df <- utils::read.delim(file = cand.file, header = T, 
                             stringsAsFactors = F, sep = "\t")


### Prepare the MAF file

# Get the same samples for the mutation file as the VEP input file
muts.df <- dplyr::filter(muts.df, candidate == TRUE &
                           vartype %in% c("NSYN", "FRS", "NFRS"))

# Merge the VEP data frame with the mutation information data frame
merged.df <- base::merge(muts.df, vep.df, by.x=c("ID", "AltAllele"), 
                   by.y=c("X.Uploaded_variation", "Allele"))

# Add the type of the mutation
merged.df <- base::transform(merged.df, Variant_Type = 
                               dplyr::case_when(
                                 RefAllele == "-" ~ "INS",
                                 AltAllele == "-" ~ "DEL",
                                 TRUE ~ "SNP"
                                 )
                             )

# Add the translational effect of variant allele
add_variant_classification <- function(conseq, type) {
  
  dplyr::case_when(
    startsWith(as.character(conseq), "synonymous_variant") == TRUE ~ "Silent", 
    startsWith(as.character(conseq), "missense_variant") == TRUE ~ "Missense_Mutation",
    startsWith(as.character(conseq), "frameshift_variant") == TRUE & 
      type == "DEL" ~ "Frame_Shift_Del",
    startsWith(as.character(conseq), "frameshift_variant") == TRUE & 
      type == "INS" ~ "Frame_Shift_Ins",
    startsWith(as.character(conseq), "protein_altering_variant") == TRUE & 
      type == "DEL"  ~ "In_Frame_Del",
    startsWith(as.character(conseq), "protein_altering_variant") == TRUE & 
      type == "INS" ~ "In_Frame_Ins",
    conseq == "inframe_deletion" ~ "In_Frame_Del",
    conseq == "inframe_insertion" ~ "In_Frame_Ins",
    startsWith(as.character(conseq), "splice_region_variant") == TRUE ~ "Splice_Region",
    startsWith(as.character(conseq), "splice_acceptor") == TRUE | 
      startsWith(as.character(conseq), "splice_donor") == TRUE ~ "Splice_Site",
    startsWith(as.character(conseq), "stop_gained") == TRUE ~ "Nonsense_Mutation",
    startsWith(as.character(conseq), "stop_lost") == TRUE ~ "Nonstop_Mutation",
    startsWith(as.character(conseq), "start_lost") == TRUE ~ "Translation_Start_Site",
    conseq == "5_prime_UTR_variant" ~ "5'UTR",
    conseq == "3_prime_UTR_variant" ~ "3'UTR",
    conseq == "up_gene_variant" ~ "5'Flank",
    conseq == "downstream_gene_variant" ~ "3'Flank",
    conseq == "intron_variant" ~ "Intron",
    startsWith(as.character(conseq), "non_coding_transcript") == TRUE ~ "RNA",  
    TRUE ~"IGR"
  )
}

merged.df <- base::transform(merged.df, Variant_Classification = 
                               add_variant_classification(
                                 conseq = merged.df$Consequence,
                                 type = merged.df$Variant_Type
                                 )
                             )

# Select the correct amino acid change
first_word <- function(AA.changes) {
  unlist(base::strsplit(AA.changes, ":"))[2]
}
merged.df$AAchange <- S4Vectors::sapply(merged.df$HGVSp, first_word)

# Change the amino acid three letter code into one letter code for readability
three_to_one_letter_code <- function(aa.change) {
  # Return immediately the change when there is no amino acid change known
  if (is.na(aa.change) == TRUE) {
    return(aa.change)
  }
  # Set the pattern for the three letter code of an amino acid
  three.code.pattern <- "[A-Z]{1}[^0-9fA-Z]{2}"
  # When there is no three letter code match, 'gregexpr' returns -1
  aa.idx <- BiocGenerics::unlist(base::gregexpr(three.code.pattern, aa.change))
  
  # Keep changing all amino acids
  while (aa.idx[1] != -1) {
    # Get the first amino acid name for replacement
    three.code <- base::substring(text = aa.change, 
                                  first = aa.idx[1], last = (aa.idx[1] + 2))
    
    # When matching the termination codon it will be replaced by an asterisk
    if (three.code == "Ter") {
      aa.change <- base::gsub(pattern = three.code, 
                              replacement = "*", x = aa.change)
    } else {
      # Use the 'a' function of seqinr for the letter code conversion
      aa.change <- base::gsub(pattern = three.code, 
                              replacement = seqinr::a(three.code), 
                              x = aa.change)
    }
    # The index will be changed by the conversion of three to one letter code
    aa.idx <- BiocGenerics::unlist(base::gregexpr(three.code.pattern, aa.change))
  }
  return(aa.change)
} 
merged.df$HGVSp_short <- S4Vectors::sapply(merged.df$AAchange, 
                                           three_to_one_letter_code)

# Rename fields to meet the MAF naming convention
merged.df <- 
  plyr::rename(merged.df, c(SYMBOL = "Hugo_Symbol", 
                            StartPos = "Start_Position", 
                            EndPos = "End_Position", 
                            RefAllele = "Reference_Allele", 
                            AltAllele = "Tumor_Seq_Allele2", 
                            caseID = "Tumor_Sample_Barcode",
                            vep.allele = "Allele"))

# Select columns needed for MAF
maf.df <-
  dplyr::select(merged.df, c(ID, Hugo_Symbol, Chromosome, Start_Position, 
                           End_Position, Reference_Allele, Tumor_Seq_Allele2, 
                           Consequence, Variant_Classification, Variant_Type, 
                           vartype, Tumor_Sample_Barcode, AAchange, 
                           HGVSp_short, ct, dataset))
