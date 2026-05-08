
library(data.table)
library(openxlsx)

############################################################
# 1) Load data
############################################################
in_file <- "/data_storage/TCR_Paired_chains.xlsx"

DT <- openxlsx::read.xlsx(in_file)
DT <- as.data.table(DT)  

# Key: CDR3alpha_CDR3beta
DT[, Merge_file := paste0(Contig_cdr3_TRA_1, "_", Contig_cdr3_TRB_1)]

# Clean Source
DT[, Source := toupper(trimws(Source))]

############################################################
# 2) Counts by Source (CSF/PBMC)
############################################################
DT_src <- DT[Source %in% c("CSF", "PBMC")]

source_counts <- DT_src[, .N, by = .(Merge_file, Source)]
source_counts <- dcast(source_counts, Merge_file ~ Source, value.var = "N", fill = 0)

if (!"CSF"  %in% names(source_counts)) source_counts[, CSF := 0L]
if (!"PBMC" %in% names(source_counts)) source_counts[, PBMC := 0L]

setnames(source_counts,
         old = c("CSF", "PBMC"),
         new = c("CSF_Count_Number", "PBMC_Count_Number"))

############################################################
# 3) Unique patient counts by Source (CSF/PBMC)
############################################################
patient_counts_by_source <- DT_src[, .(Unique_Patient_Count = uniqueN(NN_IshId_HLA)),
                                   by = .(Merge_file, Source)]
patient_counts_by_source <- dcast(patient_counts_by_source,
                                  Merge_file ~ Source,
                                  value.var = "Unique_Patient_Count",
                                  fill = 0)

if (!"CSF"  %in% names(patient_counts_by_source)) patient_counts_by_source[, CSF := 0L]
if (!"PBMC" %in% names(patient_counts_by_source)) patient_counts_by_source[, PBMC := 0L]

setnames(patient_counts_by_source,
         old = c("CSF", "PBMC"),
         new = c("Number of merged in Patient with CSF",
                 "Number of merged in Patient with PBMC"))

############################################################
# 4) Unique patients per diagnosis group (MS/CTRL/...)
############################################################
tmp_diag <- unique(DT[, .(Merge_file, NN_Diagnosis.group, NN_IshId_HLA)])
diag_counts <- tmp_diag[, .(Patient_Count = .N), by = .(Merge_file, NN_Diagnosis.group)]

phenotype_patient_counts <- dcast(diag_counts,
                                  Merge_file ~ NN_Diagnosis.group,
                                  value.var = "Patient_Count",
                                  fill = 0)

############################################################
# 4b) NEW: Unique MS/CTR IDs per Merge_file (Patient.IDs_TCR.study)
############################################################

patient_ids_study_by_combo <- DT[
  !is.na(Patient.IDs_TCR.study) & Patient.IDs_TCR.study != "",
  .(Patient.IDs_TCR.study = paste(sort(unique(Patient.IDs_TCR.study)), collapse = ",")),
  by = Merge_file
]

############################################################
# 5) MAIN SUMMARY: CDR3 + VDJ for BOTH CHAINS
############################################################
TCR_summary <- DT[, .(
  Count        = .N,
  Num_Patients = uniqueN(NN_IshId_HLA),
  
  # CDR3
  CDR3_ALPHA = Contig_cdr3_TRA_1[1],
  CDR3_BETA  = Contig_cdr3_TRB_1[1],
  
  # Contig gene calls
  TRAV = paste(unique(Contig_v_gene_TRA_1), collapse = ","),
  TRAJ = paste(unique(Contig_j_gene_TRA_1), collapse = ","),
  TRBV = paste(unique(Contig_v_gene_TRB_1), collapse = ","),
  TRBJ = paste(unique(Contig_j_gene_TRB_1), collapse = ","),
  
  # AIRR gene calls
  Airr_TRAV = paste(unique(Airr_v_call_TRA_1), collapse = ","),
  Airr_TRAJ = paste(unique(Airr_j_call_TRA_1), collapse = ","),
  Airr_TRBV = paste(unique(Airr_v_call_TRB_1), collapse = ","),
  Airr_TRBJ = paste(unique(Airr_j_call_TRB_1), collapse = ","),
  
  # Metadata lists
  NN_Diagnosis.group = paste(unique(NN_Diagnosis.group), collapse = ","),
  Airr_Sample        = paste(unique(Airr_Sample), collapse = ","),
  Source             = paste(unique(Source), collapse = ","),
  NN_IshId_HLA        = paste(unique(NN_IshId_HLA), collapse = ",")
), by = Merge_file]

############################################################
# 6) JOIN EVERYTHING
############################################################
TCR_summary_final <- merge(TCR_summary, phenotype_patient_counts, by = "Merge_file", all.x = TRUE)
TCR_summary_final <- merge(TCR_summary_final, source_counts, by = "Merge_file", all.x = TRUE)
TCR_summary_final <- merge(TCR_summary_final, patient_counts_by_source, by = "Merge_file", all.x = TRUE)

### CHANGED / ADDED ###
# add Patient.IDs_TCR.study summary per combo
TCR_summary_final <- merge(TCR_summary_final, patient_ids_study_by_combo, by = "Merge_file", all.x = TRUE)

# Replace missing numeric counts with 0
for (cc in c("CSF_Count_Number", "PBMC_Count_Number",
             "Number of merged in Patient with CSF",
             "Number of merged in Patient with PBMC")) {
  if (cc %in% names(TCR_summary_final)) {
    TCR_summary_final[is.na(get(cc)), (cc) := 0L]
  }
}


############################################################
# 7) Create Combo_All and move it to FIRST column
############################################################
TCR_summary_final[, Combo_All :=
                    paste(CDR3_ALPHA, CDR3_BETA, TRAV, TRAJ, TRBV, TRBJ, sep = "_")
]

############################################################
# 8) Rename key column + reorder columns
############################################################
setnames(TCR_summary_final, "Merge_file", "CDR3_ALPHA_CDR3_BETA")

# Put Combo_All as the FIRST column (keep everything else)
setcolorder(TCR_summary_final, c("Combo_All", setdiff(names(TCR_summary_final), "Combo_All")))
head(TCR_summary_final)
view(TCR_summary_final)
############################################################
# 9) Write Excel
############################################################
out_file <- "/data_storage/Summary_file/TCR_summary_CDR3alpha_CDR3beta_VDJ_bothChains_with_CSF_PBMC.xlsx"

openxlsx::write.xlsx(as.data.frame(TCR_summary_final), file = out_file)

