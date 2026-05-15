library(data.table)

harmonize_prs <- function(ref_path, aux_path, output_ref_path, output_aux_path) {

  ref <- fread(ref_path)
  aux <- fread(aux_path)
  
  aux[, orig_order := .I]
  
  merged <- merge(aux, ref[, .(SNP, ref_A1 = A1, ref_A2 = A2)], by = "SNP", sort = FALSE)
  
  ambiguous <- (merged$A1 == "A" & merged$A2 == "T") | 
    (merged$A1 == "T" & merged$A2 == "A") |
    (merged$A1 == "C" & merged$A2 == "G") | 
    (merged$A1 == "G" & merged$A2 == "C")
  merged <- merged[!ambiguous]

  comp <- c("A" = "T", "T" = "A", "C" = "G", "G" = "C")
  
  match_mask <- merged$A1 == merged$ref_A1 & merged$A2 == merged$ref_A2
  swap_mask <- merged$A1 == merged$ref_A2 & merged$A2 == merged$ref_A1
  strand_flip_mask <- comp[merged$A1] == merged$ref_A1 & comp[merged$A2] == merged$ref_A2
  strand_swap_mask <- comp[merged$A1] == merged$ref_A2 & comp[merged$A2] == merged$ref_A1
  
  needs_beta_flip <- swap_mask | strand_swap_mask
  merged[needs_beta_flip, BETA := BETA * -1]
  
  valid_snps <- match_mask | swap_mask | strand_flip_mask | strand_swap_mask
  final <- merged[valid_snps]
  
  final[, c("A1", "A2") := list(ref_A1, ref_A2)]
  
  setorder(final, orig_order)
  final[, c("ref_A1", "ref_A2", "orig_order") := NULL]

  fwrite(ref, output_ref_path, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  fwrite(final, output_aux_path, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  
  return(list(ref = ref, harmonized_aux = final))
}

ref_prs_path <- '<PATH TO REFERENCE (Target Population) PRS (.txt)>'
other_prs_path <- '<PATH TO OTHER (Auxiliary Population) PRS (.txt)>'
out_ref <- '<PATH TO harmonized reference (target population) PRS (same as before) (.txt)>'
out_other <- '<PATH TO harmonized other (auxiliary population) PRS (corrected) (.txt)>'

harmonize_prs(ref_prs_path, other_prs_path, out_ref, out_other)

