# Script to run 2-sample MR in Macronutrients
library(MendelianRandomization)
library(TwoSampleMR)
library(rio)
library(tidyverse)

#Read exposures or download locally from https://www.thessgac.org/data
# Proteins
prot_fil <- import("Diet_Protein_GWAS_MA_SSGAC_2020_MolPsych.txt")
# add phenotype
prot_fil$Phenotype <- "protein"
#filter pvalue
prot_fil_names<- prot_fil[prot_fil$Pval < 5e-06 ,]
dim(prot_fil_names)
#change names for format
names(prot_fil_names)[names(prot_fil_names) == "rsID"] <- "SNP"
names(prot_fil_names)[names(prot_fil_names) == "CHR"] <- "Gene"
names(prot_fil_names)[names(prot_fil_names) == "POS"] <- "pos"
names(prot_fil_names)[names(prot_fil_names) == "A1"] <- "effect_allele"
names(prot_fil_names)[names(prot_fil_names) == "A2"] <- "other_allele"
names(prot_fil_names)[names(prot_fil_names) == "FREQA1_HRC"] <- "eaf"
names(prot_fil_names)[names(prot_fil_names) == "N"] <- "samplesize"
names(prot_fil_names)[names(prot_fil_names) == "Pval"] <- "pval"
names(prot_fil_names)[names(prot_fil_names) == "Beta"] <- "beta"
names(prot_fil_names)[names(prot_fil_names) == "SE"] <- "se"
# format
prot_exp_dat <- format_data(prot_fil_names, type="exposure")
head(prot_exp_dat)
# clump
prot_exp_dat <- clump_data(prot_exp_dat)

# FAT
fat_fil <- import("Diet_Fat_GWAS_MA_SSGAC_2020_MolPsych.txt")
# add phenotype
fat_fil$Phenotype <- "fat"
#filter pvalue
fat_fil_names<- fat_fil[fat_fil$Pval < 5e-06,]
dim(fat_fil_names)
#change names for format
names(fat_fil_names)[names(fat_fil_names) == "rsID"] <- "SNP"
names(fat_fil_names)[names(fat_fil_names) == "CHR"] <- "Gene"
names(fat_fil_names)[names(fat_fil_names) == "POS"] <- "pos"
names(fat_fil_names)[names(fat_fil_names) == "A1"] <- "effect_allele"
names(fat_fil_names)[names(fat_fil_names) == "A2"] <- "other_allele"
names(fat_fil_names)[names(fat_fil_names) == "FREQA1_HRC"] <- "eaf"
names(fat_fil_names)[names(fat_fil_names) == "N"] <- "samplesize"
names(fat_fil_names)[names(fat_fil_names) == "Pval"] <- "pval"
names(fat_fil_names)[names(fat_fil_names) == "Beta"] <- "beta"
names(fat_fil_names)[names(fat_fil_names) == "SE"] <- "se"
# format
fat_exp_dat <- format_data(fat_fil_names, type="exposure")
head(fat_exp_dat)
# clump
fat_exp_dat <- clump_data(fat_exp_dat)

# CARBS
carb_fil <- import("Diet_Carbohydrate_GWAS_MA_SSGAC_2020_MolPsych.txt")
# add phenotype
carb_fil$Phenotype <- "carb"
#filter pvalue
carb_fil_names<- carb_fil[carb_fil$Pval < 5e-06,]
dim(carb_fil_names)
#change names for format
names(carb_fil_names)[names(carb_fil_names) == "rsID"] <- "SNP"
names(carb_fil_names)[names(carb_fil_names) == "CHR"] <- "Gene"
names(carb_fil_names)[names(carb_fil_names) == "POS"] <- "pos"
names(carb_fil_names)[names(carb_fil_names) == "A1"] <- "effect_allele"
names(carb_fil_names)[names(carb_fil_names) == "A2"] <- "other_allele"
names(carb_fil_names)[names(carb_fil_names) == "FREQA1_HRC"] <- "eaf"
names(carb_fil_names)[names(carb_fil_names) == "N"] <- "samplesize"
names(carb_fil_names)[names(carb_fil_names) == "Pval"] <- "pval"
names(carb_fil_names)[names(carb_fil_names) == "Beta"] <- "beta"
names(carb_fil_names)[names(carb_fil_names) == "SE"] <- "se"
# format
carb_exp_dat <- format_data(carb_fil_names, type="exposure")
head(carb_exp_dat)
# clump
carb_exp_dat <- clump_data(carb_exp_dat)

# SUGAR
sugar_fil <- import("Diet_Sugar_GWAS_MA_SSGAC_2020_MolPsych.txt")
# add phenotype
sugar_fil$Phenotype <- "sugar"
#filter pvalue
sugar_fil_names<- sugar_fil[sugar_fil$Pval < 5e-06,]
dim(sugar_fil_names)
#change names for format
names(sugar_fil_names)[names(sugar_fil_names) == "rsID"] <- "SNP"
names(sugar_fil_names)[names(sugar_fil_names) == "CHR"] <- "Gene"
names(sugar_fil_names)[names(sugar_fil_names) == "POS"] <- "pos"
names(sugar_fil_names)[names(sugar_fil_names) == "A1"] <- "effect_allele"
names(sugar_fil_names)[names(sugar_fil_names) == "A2"] <- "other_allele"
names(sugar_fil_names)[names(sugar_fil_names) == "FREQA1_HRC"] <- "eaf"
names(sugar_fil_names)[names(sugar_fil_names) == "N"] <- "samplesize"
names(sugar_fil_names)[names(sugar_fil_names) == "Pval"] <- "pval"
names(sugar_fil_names)[names(sugar_fil_names) == "Beta"] <- "beta"
names(sugar_fil_names)[names(sugar_fil_names) == "SE"] <- "se"
# format
sugar_exp_dat <- format_data(sugar_fil_names, type="exposure")
head(sugar_exp_dat)
# clump
sugar_exp_dat <- clump_data(sugar_exp_dat)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Check outcomes
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #ao <- available_outcomes()
ao[grepl("hdl", ao$trait), ] # ieu-b-109 (PMID 32203549)
ao[grepl("ldl", ao$trait), ] # ieu-b-110
ao[grepl("cholesterol", ao$trait), ] # met-d-Total_C (PMID: 32114887)/ ukb-d-30690_irnt
ao[grepl("triglycerides", ao$trait), ] # ieu-b-111
ao[grepl("stroke", ao$trait), ] # ebi-a-GCST006906
ao[grepl("heart disease", ao$trait), ]
ao[grepl("diabetes", ao$trait), ] #GCST007515 , or #GCST007516 or #GCST007517

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Read outcomes for Proteins
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
chd_out_dat_prot <- extract_outcome_data(
    snps = prot_exp_dat$SNP,
    outcomes = 'ieu-a-7'
)

t2dadjbmi_out_dat_prot <- extract_outcome_data(
    snps = prot_exp_dat$SNP,
    outcomes = 'ebi-a-GCST007518'
)
t2d_out_dat_prot <- extract_outcome_data(
    snps = prot_exp_dat$SNP,
    outcomes = 'ebi-a-GCST007517'
)
fg_out_dat_prot <- extract_outcome_data(
    snps = prot_exp_dat$SNP,
    outcomes = 'ieu-b-113'
)

twohg_out_dat_prot <- extract_outcome_data(
    snps = prot_exp_dat$SNP,
    outcomes = 'ebi-a-GCST000569'
)
stroke_out_dat_prot <- extract_outcome_data(
    snps = prot_exp_dat$SNP,
    outcomes = 'ebi-a-GCST006906'
)
hdl_out_dat_prot <- extract_outcome_data(
    snps = prot_exp_dat$SNP,
    outcomes = 'ieu-b-109'
)
ldl_out_dat_prot <- extract_outcome_data(
    snps = prot_exp_dat$SNP,
    outcomes = 'ieu-b-110'
)
skol_out_dat_prot <- extract_outcome_data(
    snps = prot_exp_dat$SNP,
    outcomes = 'met-d-Total_C'
)
stg_out_dat_prot <- extract_outcome_data(
    snps = prot_exp_dat$SNP,
    outcomes = 'ieu-b-111'
)

# CHD
# Harmonise
dat_harm_chd_prot<- harmonise_data(
    exposure_dat = prot_exp_dat,
    outcome_dat = chd_out_dat_prot)
#
mr_report(dat_harm_chd_prot,output_path = "YOUR_PATHWAY_DIRECTORY/",
          output_type = "pdf",author = "Hugo P", study = "Two Sample MR",path = system.file("reports", package = "TwoSampleMR"))

# create input
MRInput_chd_prot <- mr_input(bx = dat_harm_chd_prot$beta.exposure,
                             bxse = dat_harm_chd_prot$se.exposure,
                             by = dat_harm_chd_prot$beta.outcome,
                             byse = dat_harm_chd_prot$se.outcome)
#
raps_dat_harm_chd_prot <- mr(dat_harm_chd_prot,  method_list = c("mr_raps"))
#Run MR
MRAll_chd_prot <- mr_allmethods(MRInput_chd_prot, method = "all")
res_chd_prot<- mr(dat_harm_chd_prot)
#
sink("YOUR_PATHWAY_DIRECTORY/MRAll_chd_prot.txt")
MRAll_chd_prot
generate_odds_ratios(res_chd_prot)
sink()
p_chd_prot<- mr_plot(mr_allmethods(MRInput_chd_prot, method = "main"))
p_chd_prot +
    xlab("Protein intake (E%)") +
    ylab("Coronary heart disease (CHD)")+
    ggsave(file = "YOUR_PATHWAY_DIRECTORY/MRAll_chd_prot_plot.pdf")

# Sensitivity analyses
het_comb_chd_prot<-mr_heterogeneity(dat_harm_chd_prot) #Heterogeneity Q and Q p-val
plt_comb_chd_prot<-mr_pleiotropy_test(dat_harm_chd_prot)  #Horizontal pleiotropy
sin_comb_chd_prot<-mr_singlesnp(dat_harm_chd_prot, all_method=c("mr_ivw", "mr_egger_regression"))
p_chd_prot <- mr_forest_plot(sin_comb_chd_prot)
p_chd_prot[[1]]
dir_chd_prot <- directionality_test(dat_harm_chd_prot)
leave_comb_chd_prot<- mr_leaveoneout(dat_harm_chd_prot, parameters = default_parameters())
p_leave_chd_prot <- mr_leaveoneout_plot(leave_comb_chd_prot)
p_leave_chd_prot[[1]]
#
mrpresso_chd_prot<- run_mr_presso(dat_harm_chd_prot, NbDistribution = 10000, SignifThreshold = 0.05)
saveRDS(mrpresso_chd_prot, file = "YOUR_PATHWAY_DIRECTORY/mrpresso_chd_prot.rds")

# t2dadjbmi
# Harmonise
dat_harm_t2dadjbmi_prot<- harmonise_data(
    exposure_dat = prot_exp_dat,
    outcome_dat = t2dadjbmi_out_dat_prot)
#
mr_report(dat_harm_t2dadjbmi_prot,output_path = "YOUR_PATHWAY_DIRECTORY/",
          output_type = "pdf",author = "Hugo P", study = "Two Sample MR",path = system.file("reports", package = "TwoSampleMR"))

# create input
MRInput_t2dadjbmi_prot <- mr_input(bx = dat_harm_t2dadjbmi_prot$beta.exposure,
                                   bxse = dat_harm_t2dadjbmi_prot$se.exposure,
                                   by = dat_harm_t2dadjbmi_prot$beta.outcome,
                                   byse = dat_harm_t2dadjbmi_prot$se.outcome)
#Run MR
MRAll_t2dadjbmi_prot <- mr_allmethods(MRInput_t2dadjbmi_prot, method = "all")
res_t2dadjbmi_prot<- mr(dat_harm_t2dadjbmi_prot)
sink("YOUR_PATHWAY_DIRECTORY/MRAll_t2dadjbmi_prot.txt")
MRAll_t2dadjbmi_prot
generate_odds_ratios(res_t2dadjbmi_prot)
sink()
p_t2dadjbmi_prot<- mr_plot(mr_allmethods(MRInput_t2dadjbmi_prot, method = "main"))
p_t2dadjbmi_prot +
    xlab("Protein intake (E%)") +
    ylab("T2D adjusted for BMI (T2DadjBMI)") +
    ggsave(file = "YOUR_PATHWAY_DIRECTORY/MRAll_t2dadjbmi_prot_plot.pdf")

# Sensitivity analyses
het_comb_t2dadjbmi_prot<-mr_heterogeneity(dat_harm_t2dadjbmi_prot) #Heterogeneity Q and Q p-val
plt_comb_t2dadjbmi_prot<-mr_pleiotropy_test(dat_harm_t2dadjbmi_prot)  #Horizontal pleiotropy
sin_comb_t2dadjbmi_prot<-mr_singlesnp(dat_harm_t2dadjbmi_prot, all_method=c("mr_ivw", "mr_egger_regression"))
p_t2dadjbmi_prot <- mr_forest_plot(sin_comb_t2dadjbmi_prot)
p_t2dadjbmi_prot[[1]]
dir_t2dadjbmi_prot <- directionality_test(dat_harm_t2dadjbmi_prot)
leave_comb_t2dadjbmi_prot<- mr_leaveoneout(dat_harm_t2dadjbmi_prot, parameters = default_parameters())
p_leave_t2dadjbmi_prot <- mr_leaveoneout_plot(leave_comb_t2dadjbmi_prot)
p_leave_t2dadjbmi_prot[[1]]
#
mrpresso_t2dadjbmi_prot<- run_mr_presso(dat_harm_t2dadjbmi_prot, NbDistribution = 10000, SignifThreshold = 0.05)
saveRDS(mrpresso_t2dadjbmi_prot, file = "YOUR_PATHWAY_DIRECTORY/mrpresso_t2dadjbmi_prot.rds")
#  Not enough intrumental variables

# T2D
# Harmonise
dat_harm_t2d_prot<- harmonise_data(
    exposure_dat = prot_exp_dat,
    outcome_dat = t2d_out_dat_prot)
#
mr_report(dat_harm_t2d_prot,output_path = "YOUR_PATHWAY_DIRECTORY/",
          output_type = "pdf",author = "Hugo P", study = "Two Sample MR",path = system.file("reports", package = "TwoSampleMR"))

# create input
MRInput_t2d_prot <- mr_input(bx = dat_harm_t2d_prot$beta.exposure,
                             bxse = dat_harm_t2d_prot$se.exposure,
                             by = dat_harm_t2d_prot$beta.outcome,
                             byse = dat_harm_t2d_prot$se.outcome)
#Run MR
MRAll_t2d_prot <- mr_allmethods(MRInput_t2d_prot, method = "all")
res_t2d_prot<- mr(dat_harm_t2d_prot)
sink("YOUR_PATHWAY_DIRECTORY/MRAll_t2d_prot.txt")
MRAll_t2d_prot
generate_odds_ratios(res_t2d_prot)
sink()
p_t2d_prot<- mr_plot(mr_allmethods(MRInput_t2d_prot, method = "main"))
p_t2d_prot +
    xlab("Protein intake (E%)") +
    ylab("T2D adjusted for BMI (T2D)") +
    ggsave(file = "YOUR_PATHWAY_DIRECTORY/MRAll_t2d_prot_plot.pdf")

# Sensitivity analyses
het_comb_t2d_prot<-mr_heterogeneity(dat_harm_t2d_prot) #Heterogeneity Q and Q p-val
plt_comb_t2d_prot<-mr_pleiotropy_test(dat_harm_t2d_prot)  #Horizontal pleiotropy
sin_comb_t2d_prot<-mr_singlesnp(dat_harm_t2d_prot, all_method=c("mr_ivw", "mr_egger_regression"))
p_t2d_prot <- mr_forest_plot(sin_comb_t2d_prot)
p_t2d_prot[[1]]
dir_t2d_prot <- directionality_test(dat_harm_t2d_prot)
leave_comb_t2d_prot<- mr_leaveoneout(dat_harm_t2d_prot, parameters = default_parameters())
p_leave_t2d_prot <- mr_leaveoneout_plot(leave_comb_t2d_prot)
p_leave_t2d_prot[[1]]
#
mrpresso_t2d_prot<- run_mr_presso(dat_harm_t2d_prot, NbDistribution = 10000, SignifThreshold = 0.05)
saveRDS(mrpresso_t2d_prot, file = "YOUR_PATHWAY_DIRECTORY/mrpresso_t2d_prot.rds")

# Fasting glucose
# Harmonise
dat_harm_fg_prot<- harmonise_data(
    exposure_dat = prot_exp_dat,
    outcome_dat = fg_out_dat_prot)
#
mr_report(dat_harm_fg_prot,output_path = "YOUR_PATHWAY_DIRECTORY/",
          output_type = "pdf",author = "Hugo P", study = "Two Sample MR",path = system.file("reports", package = "TwoSampleMR"))

# create input
MRInput_fg_prot <- mr_input(bx = dat_harm_fg_prot$beta.exposure,
                            bxse = dat_harm_fg_prot$se.exposure,
                            by = dat_harm_fg_prot$beta.outcome,
                            byse = dat_harm_fg_prot$se.outcome)
#Run MR
MRAll_fg_prot <- mr_allmethods(MRInput_fg_prot, method = "all")
sink("YOUR_PATHWAY_DIRECTORY/MRAll_fg_prot.txt")
MRAll_fg_prot
sink()
p_fg_prot<- mr_plot(mr_allmethods(MRInput_fg_prot, method = "main"))
p_fg_prot +
    xlab("Protein intake (E%)") +
    ylab("Fasting glucose (FG)")+
    ggsave(file = "YOUR_PATHWAY_DIRECTORY/MRAll_fg_prot_plot.pdf")

# Sensitivity analyses
het_comb_fg_prot<-mr_heterogeneity(dat_harm_fg_prot) #Heterogeneity Q and Q p-val
plt_comb_fg_prot<-mr_pleiotropy_test(dat_harm_fg_prot)  #Horizontal pleiotropy
sin_comb_fg_prot<-mr_singlesnp(dat_harm_fg_prot, all_method=c("mr_ivw", "mr_egger_regression"))
p_fg_prot <- mr_forest_plot(sin_comb_fg_prot)
p_fg_prot[[1]]
dir_fg_prot <- directionality_test(dat_harm_fg_prot)
leave_comb_fg_prot<- mr_leaveoneout(dat_harm_fg_prot, parameters = default_parameters())
p_leave_fg_prot <- mr_leaveoneout_plot(leave_comb_fg_prot)
p_leave_fg_prot[[1]]
#
mrpresso_fg_prot<- run_mr_presso(dat_harm_fg_prot, NbDistribution = 10000, SignifThreshold = 0.05)
saveRDS(mrpresso_fg_prot, file = "YOUR_PATHWAY_DIRECTORY/mrpresso_fg_prot.rds")
#  Not enough instrumental variables

# Two Hour glucose
# Harmonise
dat_harm_twohg_prot<- harmonise_data(
    exposure_dat = prot_exp_dat,
    outcome_dat = twohg_out_dat_prot)
#
mr_report(dat_harm_twohg_prot,output_path = "YOUR_PATHWAY_DIRECTORY/",
          output_type = "pdf",author = "Hugo P", study = "Two Sample MR",path = system.file("reports", package = "TwoSampleMR"))

# create input
MRInput_twohg_prot <- mr_input(bx = dat_harm_twohg_prot$beta.exposure,
                               bxse = dat_harm_twohg_prot$se.exposure,
                               by = dat_harm_twohg_prot$beta.outcome,
                               byse = dat_harm_twohg_prot$se.outcome)
#Run MR
MRAll_twohg_prot <- mr_allmethods(MRInput_twohg_prot, method = "all")
sink("YOUR_PATHWAY_DIRECTORY/MRAll_twohg_prot.txt")
MRAll_twohg_prot
sink()
p_twohg_prot<- mr_plot(mr_allmethods(MRInput_twohg_prot, method = "main"))
p_twohg_prot +
    xlab("Protein intake (E%)") +
    ylab("Two-hour glucose (2hG)")+
    ggsave(file = "YOUR_PATHWAY_DIRECTORY/MRAll_twohg_prot_plot.pdf")

# Sensitivity analyses
het_comb_twohg_prot<-mr_heterogeneity(dat_harm_twohg_prot) #Heterogeneity Q and Q p-val
plt_comb_twohg_prot<-mr_pleiotropy_test(dat_harm_twohg_prot)  #Horizontal pleiotropy
sin_comb_twohg_prot<-mr_singlesnp(dat_harm_twohg_prot, all_method=c("mr_ivw", "mr_egger_regression"))
p_twohg_prot <- mr_forest_plot(sin_comb_twohg_prot)
p_twohg_prot[[1]]
dir_twohg_prot <- directionality_test(dat_harm_twohg_prot)
leave_comb_twohg_prot<- mr_leaveoneout(dat_harm_twohg_prot, parameters = default_parameters())
p_leave_twohg_prot <- mr_leaveoneout_plot(leave_comb_twohg_prot)
p_leave_twohg_prot[[1]]
#
mrpresso_twohg_prot<- run_mr_presso(dat_harm_twohg_prot, NbDistribution = 10000, SignifThreshold = 0.05)
saveRDS(mrpresso_twohg_prot, file = "YOUR_PATHWAY_DIRECTORY/mrpresso_twohg_prot.rds")
#  Not enough instrumental variables

# Stroke
# Harmonise
dat_harm_stroke_prot<- harmonise_data(
    exposure_dat = prot_exp_dat,
    outcome_dat = stroke_out_dat_prot)
#
mr_report(dat_harm_stroke_prot,output_path = "YOUR_PATHWAY_DIRECTORY/",
          output_type = "pdf",author = "Hugo P", study = "Two Sample MR",path = system.file("reports", package = "TwoSampleMR"))

# create input
MRInput_stroke_prot <- mr_input(bx = dat_harm_stroke_prot$beta.exposure,
                                bxse = dat_harm_stroke_prot$se.exposure,
                                by = dat_harm_stroke_prot$beta.outcome,
                                byse = dat_harm_stroke_prot$se.outcome)
#Run MR
MRAll_stroke_prot <- mr_allmethods(MRInput_stroke_prot, method = "all")
res_stroke_prot<- mr(dat_harm_stroke_prot)
sink("YOUR_PATHWAY_DIRECTORY/MRAll_stroke_prot.txt")
MRAll_stroke_prot
generate_odds_ratios(res_stroke_prot)
sink()
p_stroke_prot<- mr_plot(mr_allmethods(MRInput_stroke_prot, method = "main"))
p_stroke_prot +
    xlab("Protein intake (E%)") +
    ylab("Stroke")+
    ggsave(file = "YOUR_PATHWAY_DIRECTORY/MRAll_stroke_prot_plot.pdf")

# Sensitivity analyses
het_comb_stroke_prot<-mr_heterogeneity(dat_harm_stroke_prot) #Heterogeneity Q and Q p-val
plt_comb_stroke_prot<-mr_pleiotropy_test(dat_harm_stroke_prot)  #Horizontal pleiotropy
sin_comb_stroke_prot<-mr_singlesnp(dat_harm_stroke_prot, all_method=c("mr_ivw", "mr_egger_regression"))
p_stroke_prot <- mr_forest_plot(sin_comb_stroke_prot)
p_stroke_prot[[1]]
dir_stroke_prot <- directionality_test(dat_harm_stroke_prot)
leave_comb_stroke_prot<- mr_leaveoneout(dat_harm_stroke_prot, parameters = default_parameters())
p_leave_stroke_prot <- mr_leaveoneout_plot(leave_comb_stroke_prot)
p_leave_stroke_prot[[1]]
#
mrpresso_stroke_prot<- run_mr_presso(dat_harm_stroke_prot, NbDistribution = 10000, SignifThreshold = 0.05)
saveRDS(mrpresso_stroke_prot, file = "YOUR_PATHWAY_DIRECTORY/mrpresso_stroke_prot.rds")

# HDL
# Harmonise
dat_harm_hdl_prot<- harmonise_data(
    exposure_dat = prot_exp_dat,
    outcome_dat = hdl_out_dat_prot)
#
mr_report(dat_harm_hdl_prot,output_path = "YOUR_PATHWAY_DIRECTORY/",
          output_type = "pdf",author = "Hugo P", study = "Two Sample MR",path = system.file("reports", package = "TwoSampleMR"))

# create input
MRInput_hdl_prot <- mr_input(bx = dat_harm_hdl_prot$beta.exposure,
                             bxse = dat_harm_hdl_prot$se.exposure,
                             by = dat_harm_hdl_prot$beta.outcome,
                             byse = dat_harm_hdl_prot$se.outcome)
#Run MR
MRAll_hdl_prot <- mr_allmethods(MRInput_hdl_prot, method = "all")
sink("YOUR_PATHWAY_DIRECTORY/MRAll_hdl_prot.txt")
MRAll_hdl_prot
sink()
p_hdl_prot<- mr_plot(mr_allmethods(MRInput_hdl_prot, method = "main"))
p_hdl_prot +
    xlab("Protein intake (E%)") +
    ylab("hdl")+
    ggsave(file = "YOUR_PATHWAY_DIRECTORY/MRAll_hdl_prot_plot.pdf")

# Sensitivity analyses
het_comb_hdl_prot<-mr_heterogeneity(dat_harm_hdl_prot) #Heterogeneity Q and Q p-val
plt_comb_hdl_prot<-mr_pleiotropy_test(dat_harm_hdl_prot)  #Horizontal pleiotropy
sin_comb_hdl_prot<-mr_singlesnp(dat_harm_hdl_prot, all_method=c("mr_ivw", "mr_egger_regression"))
p_hdl_prot <- mr_forest_plot(sin_comb_hdl_prot)
p_hdl_prot[[1]]
dir_hdl_prot <- directionality_test(dat_harm_hdl_prot)
leave_comb_hdl_prot<- mr_leaveoneout(dat_harm_hdl_prot, parameters = default_parameters())
p_leave_hdl_prot <- mr_leaveoneout_plot(leave_comb_hdl_prot)
p_leave_hdl_prot[[1]]
#
mrpresso_hdl_prot<- run_mr_presso(dat_harm_hdl_prot, NbDistribution = 10000, SignifThreshold = 0.05)
saveRDS(mrpresso_hdl_prot, file = "YOUR_PATHWAY_DIRECTORY/mrpresso_hdl_prot.rds")

# LDL
# Harmonise
dat_harm_ldl_prot<- harmonise_data(
    exposure_dat = prot_exp_dat,
    outcome_dat = ldl_out_dat_prot)
#
mr_report(dat_harm_ldl_prot,output_path = "YOUR_PATHWAY_DIRECTORY/",
          output_type = "pdf",author = "Hugo P", study = "Two Sample MR",path = system.file("reports", package = "TwoSampleMR"))

# create input
MRInput_ldl_prot <- mr_input(bx = dat_harm_ldl_prot$beta.exposure,
                             bxse = dat_harm_ldl_prot$se.exposure,
                             by = dat_harm_ldl_prot$beta.outcome,
                             byse = dat_harm_ldl_prot$se.outcome)
#Run MR
MRAll_ldl_prot <- mr_allmethods(MRInput_ldl_prot, method = "all")
sink("YOUR_PATHWAY_DIRECTORY/MRAll_ldl_prot.txt")
MRAll_ldl_prot
sink()
p_ldl_prot<- mr_plot(mr_allmethods(MRInput_ldl_prot, method = "main"))
p_ldl_prot +
    xlab("Protein intake (E%)") +
    ylab("ldl")+
    ggsave(file = "YOUR_PATHWAY_DIRECTORY/MRAll_ldl_prot_plot.pdf")

# Sensitivity analyses
het_comb_ldl_prot<-mr_heterogeneity(dat_harm_ldl_prot) #Heterogeneity Q and Q p-val
plt_comb_ldl_prot<-mr_pleiotropy_test(dat_harm_ldl_prot)  #Horizontal pleiotropy
sin_comb_ldl_prot<-mr_singlesnp(dat_harm_ldl_prot, all_method=c("mr_ivw", "mr_egger_regression"))
p_ldl_prot <- mr_forest_plot(sin_comb_ldl_prot)
p_ldl_prot[[1]]
dir_ldl_prot <- directionality_test(dat_harm_ldl_prot)
leave_comb_ldl_prot<- mr_leaveoneout(dat_harm_ldl_prot, parameters = default_parameters())
p_leave_ldl_prot <- mr_leaveoneout_plot(leave_comb_ldl_prot)
p_leave_ldl_prot[[1]]
#
mrpresso_ldl_prot<- run_mr_presso(dat_harm_ldl_prot, NbDistribution = 10000, SignifThreshold = 0.05)
saveRDS(mrpresso_ldl_prot, file = "YOUR_PATHWAY_DIRECTORY/mrpresso_ldl_prot.rds")

# CHOLESTEROL
# Harmonise
dat_harm_skol_prot<- harmonise_data(
    exposure_dat = prot_exp_dat,
    outcome_dat = skol_out_dat_prot)
#
mr_report(dat_harm_skol_prot,output_path = "YOUR_PATHWAY_DIRECTORY/",
          output_type = "pdf",author = "Hugo P", study = "Two Sample MR",path = system.file("reports", package = "TwoSampleMR"))

# create input
MRInput_skol_prot <- mr_input(bx = dat_harm_skol_prot$beta.exposure,
                              bxse = dat_harm_skol_prot$se.exposure,
                              by = dat_harm_skol_prot$beta.outcome,
                              byse = dat_harm_skol_prot$se.outcome)
#Run MR
MRAll_skol_prot <- mr_allmethods(MRInput_skol_prot, method = "all")
sink("YOUR_PATHWAY_DIRECTORY/MRAll_skol_prot.txt")
MRAll_skol_prot
sink()
p_skol_prot<- mr_plot(mr_allmethods(MRInput_skol_prot, method = "main"))
p_skol_prot +
    xlab("Protein intake (E%)") +
    ylab("skol")+
    ggsave(file = "YOUR_PATHWAY_DIRECTORY/MRAll_skol_prot_plot.pdf")

# Sensitivity analyses
het_comb_skol_prot<-mr_heterogeneity(dat_harm_skol_prot) #Heterogeneity Q and Q p-val
plt_comb_skol_prot<-mr_pleiotropy_test(dat_harm_skol_prot)  #Horizontal pleiotropy
sin_comb_skol_prot<-mr_singlesnp(dat_harm_skol_prot, all_method=c("mr_ivw", "mr_egger_regression"))
p_skol_prot <- mr_forest_plot(sin_comb_skol_prot)
p_skol_prot[[1]]
dir_skol_prot <- directionality_test(dat_harm_skol_prot)
leave_comb_skol_prot<- mr_leaveoneout(dat_harm_skol_prot, parameters = default_parameters())
p_leave_skol_prot <- mr_leaveoneout_plot(leave_comb_skol_prot)
p_leave_skol_prot[[1]]
#
mrpresso_skol_prot<- run_mr_presso(dat_harm_skol_prot, NbDistribution = 10000, SignifThreshold = 0.05)
saveRDS(mrpresso_skol_prot, file = "YOUR_PATHWAY_DIRECTORY/mrpresso_skol_prot.rds")

# TRIGLYCERIDES
# Harmonise
dat_harm_stg_prot<- harmonise_data(
    exposure_dat = prot_exp_dat,
    outcome_dat = stg_out_dat_prot)
#
mr_report(dat_harm_stg_prot,output_path = "YOUR_PATHWAY_DIRECTORY/",
          output_type = "pdf",author = "Hugo P", study = "Two Sample MR",path = system.file("reports", package = "TwoSampleMR"))

# create input
MRInput_stg_prot <- mr_input(bx = dat_harm_stg_prot$beta.exposure,
                             bxse = dat_harm_stg_prot$se.exposure,
                             by = dat_harm_stg_prot$beta.outcome,
                             byse = dat_harm_stg_prot$se.outcome)
#Run MR
MRAll_stg_prot <- mr_allmethods(MRInput_stg_prot, method = "all")
sink("YOUR_PATHWAY_DIRECTORY/MRAll_stg_prot.txt")
MRAll_stg_prot
sink()
p_stg_prot<- mr_plot(mr_allmethods(MRInput_stg_prot, method = "main"))
p_stg_prot +
    xlab("Protein intake (E%)") +
    ylab("stg")+
    ggsave(file = "YOUR_PATHWAY_DIRECTORY/MRAll_stg_prot_plot.pdf")

# Sensitivity analyses
het_comb_stg_prot<-mr_heterogeneity(dat_harm_stg_prot) #Heterogeneity Q and Q p-val
plt_comb_stg_prot<-mr_pleiotropy_test(dat_harm_stg_prot)  #Horizontal pleiotropy
sin_comb_stg_prot<-mr_singlesnp(dat_harm_stg_prot, all_method=c("mr_ivw", "mr_egger_regression"))
p_stg_prot <- mr_forest_plot(sin_comb_stg_prot)
p_stg_prot[[1]]
dir_stg_prot <- directionality_test(dat_harm_stg_prot)
leave_comb_stg_prot<- mr_leaveoneout(dat_harm_stg_prot, parameters = default_parameters())
p_leave_stg_prot <- mr_leaveoneout_plot(leave_comb_stg_prot)
p_leave_stg_prot[[1]]
#
mrpresso_stg_prot<- run_mr_presso(dat_harm_stg_prot, NbDistribution = 10000, SignifThreshold = 0.05)
saveRDS(mrpresso_stg_prot, file = "YOUR_PATHWAY_DIRECTORY/mrpresso_stg_prot.rds")

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#Read outcomes for Fat
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
chd_out_dat_fat <- extract_outcome_data(
    snps = fat_exp_dat$SNP,
    outcomes = 'ieu-a-7'
)

t2dadjbmi_out_dat_fat <- extract_outcome_data(
    snps = fat_exp_dat$SNP,
    outcomes = 'ebi-a-GCST007518'
)
t2d_out_dat_fat <- extract_outcome_data(
    snps = fat_exp_dat$SNP,
    outcomes = 'ebi-a-GCST007517'
)
fg_out_dat_fat <- extract_outcome_data(
    snps = fat_exp_dat$SNP,
    outcomes = 'ieu-b-113'
)
twohg_out_dat_fat <- extract_outcome_data(
    snps = fat_exp_dat$SNP,
    outcomes = 'ebi-a-GCST000569'
)
stroke_out_dat_fat <- extract_outcome_data(
    snps = prot_exp_dat$SNP,
    outcomes = 'ebi-a-GCST006906'
)
hdl_out_dat_fat <- extract_outcome_data(
    snps = fat_exp_dat$SNP,
    outcomes = 'ieu-b-109'
)
ldl_out_dat_fat <- extract_outcome_data(
    snps = fat_exp_dat$SNP,
    outcomes = 'ieu-b-110'
)
skol_out_dat_fat <- extract_outcome_data(
    snps = fat_exp_dat$SNP,
    outcomes = 'met-d-Total_C'
)
stg_out_dat_fat <- extract_outcome_data(
    snps = fat_exp_dat$SNP,
    outcomes = 'ieu-b-111'
)

# CHD
# Harmonise
dat_harm_chd_fat<- harmonise_data(
    exposure_dat = fat_exp_dat,
    outcome_dat = chd_out_dat_fat)
#
mr_report(dat_harm_chd_fat,output_path = "YOUR_PATHWAY_DIRECTORY/",
          output_type = "pdf",author = "Hugo P", study = "Two Sample MR",path = system.file("reports", package = "TwoSampleMR"))

# create input
MRInput_chd_fat <- mr_input(bx = dat_harm_chd_fat$beta.exposure,
                            bxse = dat_harm_chd_fat$se.exposure,
                            by = dat_harm_chd_fat$beta.outcome,
                            byse = dat_harm_chd_fat$se.outcome)
#Run MR
MRAll_chd_fat <- mr_allmethods(MRInput_chd_fat, method = "all")
res_chd_fat<- mr(dat_harm_chd_fat)
sink("YOUR_PATHWAY_DIRECTORY/MRAll_chd_fat.txt")
MRAll_chd_fat
generate_odds_ratios(res_chd_fat)
sink()
p_chd_fat<- mr_plot(mr_allmethods(MRInput_chd_fat, method = "main"))
p_chd_fat +
    xlab("Fat intake (E%)") +
    ylab("Coronary heart disease (CHD)")+
    ggsave(file = "YOUR_PATHWAY_DIRECTORY/MRAll_chd_fat_plot.pdf")

# Sensitivity analyses
het_comb_chd_fat<-mr_heterogeneity(dat_harm_chd_fat) #Heterogeneity Q and Q p-val
plt_comb_chd_fat<-mr_pleiotropy_test(dat_harm_chd_fat)  #Horizontal pleiotropy
sin_comb_chd_fat<-mr_singlesnp(dat_harm_chd_fat, all_method=c("mr_ivw", "mr_egger_regression"))
p_chd_fat <- mr_forest_plot(sin_comb_chd_fat)
p_chd_fat[[1]]
dir_chd_fat <- directionality_test(dat_harm_chd_fat)
leave_comb_chd_fat<- mr_leaveoneout(dat_harm_chd_fat, parameters = default_parameters())
p_leave_chd_fat <- mr_leaveoneout_plot(leave_comb_chd_fat)
p_leave_chd_fat[[1]]
#
mrpresso_chd_fat<- run_mr_presso(dat_harm_chd_fat, NbDistribution = 10000, SignifThreshold = 0.05)
saveRDS(mrpresso_chd_fat, file = "YOUR_PATHWAY_DIRECTORY/mrpresso_chd_fat.rds")

# t2dadjbmi
# Harmonise
dat_harm_t2dadjbmi_fat<- harmonise_data(
    exposure_dat = fat_exp_dat,
    outcome_dat = t2dadjbmi_out_dat_fat)
#
mr_report(dat_harm_t2dadjbmi_fat,output_path = "YOUR_PATHWAY_DIRECTORY/",
          output_type = "pdf",author = "Hugo P", study = "Two Sample MR",path = system.file("reports", package = "TwoSampleMR"))

# create input
MRInput_t2dadjbmi_fat <- mr_input(bx = dat_harm_t2dadjbmi_fat$beta.exposure,
                                  bxse = dat_harm_t2dadjbmi_fat$se.exposure,
                                  by = dat_harm_t2dadjbmi_fat$beta.outcome,
                                  byse = dat_harm_t2dadjbmi_fat$se.outcome)
#Run MR
MRAll_t2dadjbmi_fat <- mr_allmethods(MRInput_t2dadjbmi_fat, method = "all")
res_t2dadjbmi_fat<- mr(dat_harm_t2dadjbmi_fat)
sink("YOUR_PATHWAY_DIRECTORY/MRAll_t2dadjbmi_fat.txt")
MRAll_t2dadjbmi_fat
generate_odds_ratios(res_t2dadjbmi_fat)
sink()
p_t2dadjbmi_fat<- mr_plot(mr_allmethods(MRInput_t2dadjbmi_fat, method = "main"))
p_t2dadjbmi_fat +
    xlab("Fat intake (E%)") +
    ylab("T2D adjusted for BMI (T2DadjBMI)")+
    ggsave(file = "YOUR_PATHWAY_DIRECTORY/MRAll_t2dadjbmi_fat_plot.pdf")

# Sensitivity analyses
het_comb_t2dadjbmi_fat<-mr_heterogeneity(dat_harm_t2dadjbmi_fat) #Heterogeneity Q and Q p-val
plt_comb_t2dadjbmi_fat<-mr_pleiotropy_test(dat_harm_t2dadjbmi_fat)  #Horizontal pleiotropy
sin_comb_t2dadjbmi_fat<-mr_singlesnp(dat_harm_t2dadjbmi_fat, all_method=c("mr_ivw", "mr_egger_regression"))
p_t2dadjbmi_fat <- mr_forest_plot(sin_comb_t2dadjbmi_fat)
p_t2dadjbmi_fat[[1]]
dir_t2dadjbmi_fat <- directionality_test(dat_harm_t2dadjbmi_fat)
leave_comb_t2dadjbmi_fat<- mr_leaveoneout(dat_harm_t2dadjbmi_fat, parameters = default_parameters())
p_leave_t2dadjbmi_fat <- mr_leaveoneout_plot(leave_comb_t2dadjbmi_fat)
p_leave_t2dadjbmi_fat[[1]]
#
mrpresso_t2dadjbmi_fat<- run_mr_presso(dat_harm_t2dadjbmi_fat, NbDistribution = 10000, SignifThreshold = 0.05)
saveRDS(mrpresso_t2dadjbmi_fat, file = "YOUR_PATHWAY_DIRECTORY/mrpresso_t2dadjbmi_fat.rds")
# Not enough intrumental variables

# T2D
# Harmonise
dat_harm_t2d_fat<- harmonise_data(
    exposure_dat = fat_exp_dat,
    outcome_dat = t2d_out_dat_fat)
#
mr_report(dat_harm_t2d_fat,output_path = "YOUR_PATHWAY_DIRECTORY/",
          output_type = "pdf",author = "Hugo P", study = "Two Sample MR",path = system.file("reports", package = "TwoSampleMR"))

# create input
MRInput_t2d_fat <- mr_input(bx = dat_harm_t2d_fat$beta.exposure,
                            bxse = dat_harm_t2d_fat$se.exposure,
                            by = dat_harm_t2d_fat$beta.outcome,
                            byse = dat_harm_t2d_fat$se.outcome)
#Run MR
MRAll_t2d_fat <- mr_allmethods(MRInput_t2d_fat, method = "all")
res_t2d_fat<- mr(dat_harm_t2d_fat)
sink("YOUR_PATHWAY_DIRECTORY/MRAll_t2d_fat.txt")
MRAll_t2d_fat
generate_odds_ratios(res_t2d_fat)
sink()
p_t2d_fat<- mr_plot(mr_allmethods(MRInput_t2d_fat, method = "main"))
p_t2d_fat +
    xlab("fatein intake (E%)") +
    ylab("T2D adjusted for BMI (T2D)") +
    ggsave(file = "YOUR_PATHWAY_DIRECTORY/MRAll_t2d_fat_plot.pdf")

# Sensitivity analyses
het_comb_t2d_fat<-mr_heterogeneity(dat_harm_t2d_fat) #Heterogeneity Q and Q p-val
plt_comb_t2d_fat<-mr_pleiotropy_test(dat_harm_t2d_fat)  #Horizontal pleiotropy
sin_comb_t2d_fat<-mr_singlesnp(dat_harm_t2d_fat, all_method=c("mr_ivw", "mr_egger_regression"))
p_t2d_fat <- mr_forest_plot(sin_comb_t2d_fat)
p_t2d_fat[[1]]
dir_t2d_fat <- directionality_test(dat_harm_t2d_fat)
leave_comb_t2d_fat<- mr_leaveoneout(dat_harm_t2d_fat, parameters = default_parameters())
p_leave_t2d_fat <- mr_leaveoneout_plot(leave_comb_t2d_fat)
p_leave_t2d_fat[[1]]
#
mrpresso_t2d_fat<- run_mr_presso(dat_harm_t2d_fat, NbDistribution = 10000, SignifThreshold = 0.05)
saveRDS(mrpresso_t2d_fat, file = "YOUR_PATHWAY_DIRECTORY/mrpresso_t2d_fat.rds")

# Fasting glucose
# Harmonise
dat_harm_fg_fat<- harmonise_data(
    exposure_dat = fat_exp_dat,
    outcome_dat = fg_out_dat_fat)
#
mr_report(dat_harm_fg_fat,output_path = "YOUR_PATHWAY_DIRECTORY/",
          output_type = "pdf",author = "Hugo P", study = "Two Sample MR",path = system.file("reports", package = "TwoSampleMR"))

# create input
MRInput_fg_fat <- mr_input(bx = dat_harm_fg_fat$beta.exposure,
                           bxse = dat_harm_fg_fat$se.exposure,
                           by = dat_harm_fg_fat$beta.outcome,
                           byse = dat_harm_fg_fat$se.outcome)
#Run MR
MRAll_fg_fat <- mr_allmethods(MRInput_fg_fat, method = "all")
sink("YOUR_PATHWAY_DIRECTORY/MRAll_fg_fat.txt")
MRAll_fg_fat
sink()
p_fg_fat<- mr_plot(mr_allmethods(MRInput_fg_fat, method = "main"))
p_fg_fat +
    xlab("Fat intake (E%)") +
    ylab("Fasting glucose (FG)")+
    ggsave(file = "YOUR_PATHWAY_DIRECTORY/MRAll_fg_fat_plot.pdf")

# Sensitivity analyses
het_comb_fg_fat<-mr_heterogeneity(dat_harm_fg_fat) #Heterogeneity Q and Q p-val
plt_comb_fg_fat<-mr_pleiotropy_test(dat_harm_fg_fat)  #Horizontal pleiotropy
sin_comb_fg_fat<-mr_singlesnp(dat_harm_fg_fat, all_method=c("mr_ivw", "mr_egger_regression"))
p_fg_fat <- mr_forest_plot(sin_comb_fg_fat)
p_fg_fat[[1]]
dir_fg_fat <- directionality_test(dat_harm_fg_fat)
leave_comb_fg_fat<- mr_leaveoneout(dat_harm_fg_fat, parameters = default_parameters())
p_leave_fg_fat <- mr_leaveoneout_plot(leave_comb_fg_fat)
p_leave_fg_fat[[1]]
#
mrpresso_fg_fat<- run_mr_presso(dat_harm_fg_fat, NbDistribution = 10000, SignifThreshold = 0.05)
saveRDS(mrpresso_fg_fat, file = "YOUR_PATHWAY_DIRECTORY/mrpresso_fg_fat.rds")
#  Not enough instrumental variables

# Two Hour glucose
# Harmonise
dat_harm_twohg_fat<- harmonise_data(
    exposure_dat = fat_exp_dat,
    outcome_dat = twohg_out_dat_fat)
#
mr_report(dat_harm_twohg_fat,output_path = "YOUR_PATHWAY_DIRECTORY/",
          output_type = "pdf",author = "Hugo P", study = "Two Sample MR",path = system.file("reports", package = "TwoSampleMR"))

# create input
MRInput_twohg_fat <- mr_input(bx = dat_harm_twohg_fat$beta.exposure,
                              bxse = dat_harm_twohg_fat$se.exposure,
                              by = dat_harm_twohg_fat$beta.outcome,
                              byse = dat_harm_twohg_fat$se.outcome)
#Run MR
MRAll_twohg_fat <- mr_allmethods(MRInput_twohg_fat, method = "all")
sink("YOUR_PATHWAY_DIRECTORY/MRAll_twohg_fat.txt")
MRAll_twohg_fat
sink()
p_twohg_fat<- mr_plot(mr_allmethods(MRInput_twohg_fat, method = "main"))
p_twohg_fat +
    xlab("Fat intake (E%)") +
    ylab("Two-hour glucose (2hG)")+
    ggsave(file = "YOUR_PATHWAY_DIRECTORY/MRAll_twohg_fat_plot.pdf")

# Sensitivity analyses
het_comb_twohg_fat<-mr_heterogeneity(dat_harm_twohg_fat) #Heterogeneity Q and Q p-val
plt_comb_twohg_fat<-mr_pleiotropy_test(dat_harm_twohg_fat)  #Horizontal pleiotropy
sin_comb_twohg_fat<-mr_singlesnp(dat_harm_twohg_fat, all_method=c("mr_ivw", "mr_egger_regression"))
p_twohg_fat <- mr_forest_plot(sin_comb_twohg_fat)
p_twohg_fat[[1]]
dir_twohg_fat <- directionality_test(dat_harm_twohg_fat)
leave_comb_twohg_fat<- mr_leaveoneout(dat_harm_twohg_fat, parameters = default_parameters())
p_leave_twohg_fat <- mr_leaveoneout_plot(leave_comb_twohg_fat)
p_leave_twohg_fat[[1]]
#
mrpresso_twohg_fat<- run_mr_presso(dat_harm_twohg_fat, NbDistribution = 10000, SignifThreshold = 0.05)
saveRDS(mrpresso_twohg_fat, file = "YOUR_PATHWAY_DIRECTORY/mrpresso_twohg_fat.rds")
#  Not enough instrumental variables

# Stroke
# Harmonise
dat_harm_stroke_fat<- harmonise_data(
    exposure_dat = fat_exp_dat,
    outcome_dat = stroke_out_dat_fat)
#
mr_report(dat_harm_stroke_fat,output_path = "YOUR_PATHWAY_DIRECTORY/",
          output_type = "pdf",author = "Hugo P", study = "Two Sample MR",path = system.file("reports", package = "TwoSampleMR"))

# create input
MRInput_stroke_fat <- mr_input(bx = dat_harm_stroke_fat$beta.exposure,
                               bxse = dat_harm_stroke_fat$se.exposure,
                               by = dat_harm_stroke_fat$beta.outcome,
                               byse = dat_harm_stroke_fat$se.outcome)
#Run MR
MRAll_stroke_fat <- mr_allmethods(MRInput_stroke_fat, method = "all") #Method requires data on >2 variants
res_stroke_fat<- mr(dat_harm_stroke_fat)
sink("YOUR_PATHWAY_DIRECTORY/MRAll_stroke_fat.txt")
MRAll_stroke_fat
generate_odds_ratios(res_stroke_fat)
sink()
p_stroke_fat<- mr_plot(mr_allmethods(MRInput_stroke_fat, method = "main"))
p_stroke_fat +
    xlab("Fat intake (E%)") +
    ylab("Stroke")+
    ggsave(file = "YOUR_PATHWAY_DIRECTORY/MRAll_stroke_fat_plot.pdf")

# Sensitivity analyses
het_comb_stroke_fat<-mr_heterogeneity(dat_harm_stroke_fat) #Heterogeneity Q and Q p-val
plt_comb_stroke_fat<-mr_pleiotropy_test(dat_harm_stroke_fat)  #Horizontal pleiotropy
sin_comb_stroke_fat<-mr_singlesnp(dat_harm_stroke_fat, all_method=c("mr_ivw", "mr_egger_regression"))
p_stroke_fat <- mr_forest_plot(sin_comb_stroke_fat)
p_stroke_fat[[1]]
dir_stroke_fat <- directionality_test(dat_harm_stroke_fat)
leave_comb_stroke_fat<- mr_leaveoneout(dat_harm_stroke_fat, parameters = default_parameters())
p_leave_stroke_fat <- mr_leaveoneout_plot(leave_comb_stroke_fat)
p_leave_stroke_fat[[1]]
#
mrpresso_stroke_fat<- run_mr_presso(dat_harm_stroke_fat, NbDistribution = 10000, SignifThreshold = 0.05)
saveRDS(mrpresso_stroke_fat, file = "YOUR_PATHWAY_DIRECTORY/mrpresso_stroke_fat.rds")
#  Not enough instrumental variables

# HDL
# Harmonise
dat_harm_hdl_fat<- harmonise_data(
    exposure_dat = fat_exp_dat,
    outcome_dat = hdl_out_dat_fat)
#
mr_report(dat_harm_hdl_fat,output_path = "YOUR_PATHWAY_DIRECTORY/",
          output_type = "pdf",author = "Hugo P", study = "Two Sample MR",path = system.file("reports", package = "TwoSampleMR"))

# create input
MRInput_hdl_fat <- mr_input(bx = dat_harm_hdl_fat$beta.exposure,
                            bxse = dat_harm_hdl_fat$se.exposure,
                            by = dat_harm_hdl_fat$beta.outcome,
                            byse = dat_harm_hdl_fat$se.outcome)
#Run MR
MRAll_hdl_fat <- mr_allmethods(MRInput_hdl_fat, method = "all")
res_hdl_fat<- mr(dat_harm_hdl_fat)
sink("YOUR_PATHWAY_DIRECTORY/MRAll_hdl_fat.txt")
MRAll_hdl_fat
generate_odds_ratios(res_hdl_fat)
sink()
p_hdl_fat<- mr_plot(mr_allmethods(MRInput_hdl_fat, method = "main"))
p_hdl_fat +
    xlab("fatein intake (E%)") +
    ylab("hdl")+
    ggsave(file = "YOUR_PATHWAY_DIRECTORY/MRAll_hdl_fat_plot.pdf")

# Sensitivity analyses
het_comb_hdl_fat<-mr_heterogeneity(dat_harm_hdl_fat) #Heterogeneity Q and Q p-val
plt_comb_hdl_fat<-mr_pleiotropy_test(dat_harm_hdl_fat)  #Horizontal pleiotropy
sin_comb_hdl_fat<-mr_singlesnp(dat_harm_hdl_fat, all_method=c("mr_ivw", "mr_egger_regression"))
p_hdl_fat <- mr_forest_plot(sin_comb_hdl_fat)
p_hdl_fat[[1]]
dir_hdl_fat <- directionality_test(dat_harm_hdl_fat)
leave_comb_hdl_fat<- mr_leaveoneout(dat_harm_hdl_fat, parameters = default_parameters())
p_leave_hdl_fat <- mr_leaveoneout_plot(leave_comb_hdl_fat)
p_leave_hdl_fat[[1]]
#
mrpresso_hdl_fat<- run_mr_presso(dat_harm_hdl_fat, NbDistribution = 10000, SignifThreshold = 0.05)
saveRDS(mrpresso_hdl_fat, file = "YOUR_PATHWAY_DIRECTORY/mrpresso_hdl_fat.rds")

# LDL
# Harmonise
dat_harm_ldl_fat<- harmonise_data(
    exposure_dat = fat_exp_dat,
    outcome_dat = ldl_out_dat_fat)
#
mr_report(dat_harm_ldl_fat,output_path = "YOUR_PATHWAY_DIRECTORY/",
          output_type = "pdf",author = "Hugo P", study = "Two Sample MR",path = system.file("reports", package = "TwoSampleMR"))

# create input
MRInput_ldl_fat <- mr_input(bx = dat_harm_ldl_fat$beta.exposure,
                            bxse = dat_harm_ldl_fat$se.exposure,
                            by = dat_harm_ldl_fat$beta.outcome,
                            byse = dat_harm_ldl_fat$se.outcome)
#Run MR
MRAll_ldl_fat <- mr_allmethods(MRInput_ldl_fat, method = "all")
sink("YOUR_PATHWAY_DIRECTORY/MRAll_ldl_fat.txt")
MRAll_ldl_fat
sink()
p_ldl_fat<- mr_plot(mr_allmethods(MRInput_ldl_fat, method = "main"))
p_ldl_fat +
    xlab("fatein intake (E%)") +
    ylab("ldl")+
    ggsave(file = "YOUR_PATHWAY_DIRECTORY/MRAll_ldl_fat_plot.pdf")

# Sensitivity analyses
het_comb_ldl_fat<-mr_heterogeneity(dat_harm_ldl_fat) #Heterogeneity Q and Q p-val
plt_comb_ldl_fat<-mr_pleiotropy_test(dat_harm_ldl_fat)  #Horizontal pleiotropy
sin_comb_ldl_fat<-mr_singlesnp(dat_harm_ldl_fat, all_method=c("mr_ivw", "mr_egger_regression"))
p_ldl_fat <- mr_forest_plot(sin_comb_ldl_fat)
p_ldl_fat[[1]]
dir_ldl_fat <- directionality_test(dat_harm_ldl_fat)
leave_comb_ldl_fat<- mr_leaveoneout(dat_harm_ldl_fat, parameters = default_parameters())
p_leave_ldl_fat <- mr_leaveoneout_plot(leave_comb_ldl_fat)
p_leave_ldl_fat[[1]]
#
mrpresso_ldl_fat<- run_mr_presso(dat_harm_ldl_fat, NbDistribution = 10000, SignifThreshold = 0.05)
saveRDS(mrpresso_ldl_fat, file = "YOUR_PATHWAY_DIRECTORY/mrpresso_ldl_fat.rds")

# CHOLESTEROL
# Harmonise
dat_harm_skol_fat<- harmonise_data(
    exposure_dat = fat_exp_dat,
    outcome_dat = skol_out_dat_fat)
#
mr_report(dat_harm_skol_fat,output_path = "YOUR_PATHWAY_DIRECTORY/",
          output_type = "pdf",author = "Hugo P", study = "Two Sample MR",path = system.file("reports", package = "TwoSampleMR"))

# create input
MRInput_skol_fat <- mr_input(bx = dat_harm_skol_fat$beta.exposure,
                             bxse = dat_harm_skol_fat$se.exposure,
                             by = dat_harm_skol_fat$beta.outcome,
                             byse = dat_harm_skol_fat$se.outcome)
#Run MR
MRAll_skol_fat <- mr_allmethods(MRInput_skol_fat, method = "all")
sink("YOUR_PATHWAY_DIRECTORY/MRAll_skol_fat.txt")
MRAll_skol_fat
sink()
p_skol_fat<- mr_plot(mr_allmethods(MRInput_skol_fat, method = "main"))
p_skol_fat +
    xlab("fatein intake (E%)") +
    ylab("skol")+
    ggsave(file = "YOUR_PATHWAY_DIRECTORY/MRAll_skol_fat_plot.pdf")

# Sensitivity analyses
het_comb_skol_fat<-mr_heterogeneity(dat_harm_skol_fat) #Heterogeneity Q and Q p-val
plt_comb_skol_fat<-mr_pleiotropy_test(dat_harm_skol_fat)  #Horizontal pleiotropy
sin_comb_skol_fat<-mr_singlesnp(dat_harm_skol_fat, all_method=c("mr_ivw", "mr_egger_regression"))
p_skol_fat <- mr_forest_plot(sin_comb_skol_fat)
p_skol_fat[[1]]
dir_skol_fat <- directionality_test(dat_harm_skol_fat)
leave_comb_skol_fat<- mr_leaveoneout(dat_harm_skol_fat, parameters = default_parameters())
p_leave_skol_fat <- mr_leaveoneout_plot(leave_comb_skol_fat)
p_leave_skol_fat[[1]]
#
mrpresso_skol_fat<- run_mr_presso(dat_harm_skol_fat, NbDistribution = 10000, SignifThreshold = 0.05)
saveRDS(mrpresso_skol_fat, file = "YOUR_PATHWAY_DIRECTORY/mrpresso_skol_fat.rds")

# TRIGLYCERIDES
# Harmonise
dat_harm_stg_fat<- harmonise_data(
    exposure_dat = fat_exp_dat,
    outcome_dat = stg_out_dat_fat)
#
mr_report(dat_harm_stg_fat,output_path = "YOUR_PATHWAY_DIRECTORY/",
          output_type = "pdf",author = "Hugo P", study = "Two Sample MR",path = system.file("reports", package = "TwoSampleMR"))

# create input
MRInput_stg_fat <- mr_input(bx = dat_harm_stg_fat$beta.exposure,
                            bxse = dat_harm_stg_fat$se.exposure,
                            by = dat_harm_stg_fat$beta.outcome,
                            byse = dat_harm_stg_fat$se.outcome)
#Run MR
MRAll_stg_fat <- mr_allmethods(MRInput_stg_fat, method = "all")
sink("YOUR_PATHWAY_DIRECTORY/MRAll_stg_fat.txt")
MRAll_stg_fat
sink()
p_stg_fat<- mr_plot(mr_allmethods(MRInput_stg_fat, method = "main"))
p_stg_fat +
    xlab("fatein intake (E%)") +
    ylab("stg")+
    ggsave(file = "YOUR_PATHWAY_DIRECTORY/MRAll_stg_fat_plot.pdf")

# Sensitivity analyses
het_comb_stg_fat<-mr_heterogeneity(dat_harm_stg_fat) #Heterogeneity Q and Q p-val
plt_comb_stg_fat<-mr_pleiotropy_test(dat_harm_stg_fat)  #Horizontal pleiotropy
sin_comb_stg_fat<-mr_singlesnp(dat_harm_stg_fat, all_method=c("mr_ivw", "mr_egger_regression"))
p_stg_fat <- mr_forest_plot(sin_comb_stg_fat)
p_stg_fat[[1]]
dir_stg_fat <- directionality_test(dat_harm_stg_fat)
leave_comb_stg_fat<- mr_leaveoneout(dat_harm_stg_fat, parameters = default_parameters())
p_leave_stg_fat <- mr_leaveoneout_plot(leave_comb_stg_fat)
p_leave_stg_fat[[1]]
#
mrpresso_stg_fat<- run_mr_presso(dat_harm_stg_fat, NbDistribution = 10000, SignifThreshold = 0.05)
saveRDS(mrpresso_stg_fat, file = "YOUR_PATHWAY_DIRECTORY/mrpresso_stg_fat.rds")

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#Read outcomes for Carbs
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
chd_out_dat_carb <- extract_outcome_data(
    snps = carb_exp_dat$SNP,
    outcomes = 'ieu-a-7'
)

t2dadjbmi_out_dat_carb <- extract_outcome_data(
    snps = carb_exp_dat$SNP,
    outcomes = 'ebi-a-GCST007518' #GCST007515
)
t2d_out_dat_carb <- extract_outcome_data(
    snps = carb_exp_dat$SNP,
    outcomes = 'ebi-a-GCST007517'
)
fg_out_dat_carb <- extract_outcome_data(
    snps = carb_exp_dat$SNP,
    outcomes = 'ieu-b-113'
)

twohg_out_dat_carb <- extract_outcome_data(
    snps = carb_exp_dat$SNP,
    outcomes = 'ebi-a-GCST000569'
)
stroke_out_dat_carb <- extract_outcome_data(
    snps = prot_exp_dat$SNP,
    outcomes = 'ebi-a-GCST006906'
)
hdl_out_dat_carb <- extract_outcome_data(
    snps = carb_exp_dat$SNP,
    outcomes = 'ieu-b-109'
)
ldl_out_dat_carb <- extract_outcome_data(
    snps = carb_exp_dat$SNP,
    outcomes = 'ieu-b-110'
)
skol_out_dat_carb <- extract_outcome_data(
    snps = carb_exp_dat$SNP,
    outcomes = 'met-d-Total_C'
)
stg_out_dat_carb <- extract_outcome_data(
    snps = carb_exp_dat$SNP,
    outcomes = 'ieu-b-111'
)

# CHD
# Harmonise
dat_harm_chd_carb<- harmonise_data(
    exposure_dat = carb_exp_dat,
    outcome_dat = chd_out_dat_carb)
#
mr_report(dat_harm_chd_carb,output_path = "YOUR_PATHWAY_DIRECTORY/",
          output_type = "pdf",author = "Hugo P", study = "Two Sample MR",path = system.file("reports", package = "TwoSampleMR"))

# create input
MRInput_chd_carb <- mr_input(bx = dat_harm_chd_carb$beta.exposure,
                             bxse = dat_harm_chd_carb$se.exposure,
                             by = dat_harm_chd_carb$beta.outcome,
                             byse = dat_harm_chd_carb$se.outcome)
#Run MR
MRAll_chd_carb <- mr_allmethods(MRInput_chd_carb, method = "all")
res_chd_carb<- mr(dat_harm_chd_carb)
sink("YOUR_PATHWAY_DIRECTORY/MRAll_chd_carb.txt")
MRAll_chd_carb
generate_odds_ratios(res_chd_carb)
sink()
p_chd_carb<- mr_plot(mr_allmethods(MRInput_chd_carb, method = "main"))
p_chd_carb +
    xlab("Carbohydrates intake (E%)") +
    ylab("Coronary heart disease (CHD)")+
    ggsave(file = "YOUR_PATHWAY_DIRECTORY/MRAll_chd_carb_plot.pdf")

# Sensitivity analyses
het_comb_chd_carb<-mr_heterogeneity(dat_harm_chd_carb) #Heterogeneity Q and Q p-val
plt_comb_chd_carb<-mr_pleiotropy_test(dat_harm_chd_carb)  #Horizontal pleiotropy
sin_comb_chd_carb<-mr_singlesnp(dat_harm_chd_carb, all_method=c("mr_ivw", "mr_egger_regression"))
p_chd_carb <- mr_forest_plot(sin_comb_chd_carb)
p_chd_carb[[1]]
dir_chd_carb <- directionality_test(dat_harm_chd_carb)
leave_comb_chd_carb<- mr_leaveoneout(dat_harm_chd_carb, parameters = default_parameters())
p_leave_chd_carb <- mr_leaveoneout_plot(leave_comb_chd_carb)
p_leave_chd_carb[[1]]
#
mrpresso_chd_carb<- run_mr_presso(dat_harm_chd_carb, NbDistribution = 10000, SignifThreshold = 0.05)
saveRDS(mrpresso_chd_carb, file = "YOUR_PATHWAY_DIRECTORY/mrpresso_chd_carb.rds")

# t2dadjbmi
# Harmonise
dat_harm_t2dadjbmi_carb<- harmonise_data(
    exposure_dat = carb_exp_dat,
    outcome_dat = t2dadjbmi_out_dat_carb)
#
mr_report(dat_harm_t2dadjbmi_carb,output_path = "YOUR_PATHWAY_DIRECTORY/",
          output_type = "pdf",author = "Hugo P", study = "Two Sample MR",path = system.file("reports", package = "TwoSampleMR"))

# create input
MRInput_t2dadjbmi_carb <- mr_input(bx = dat_harm_t2dadjbmi_carb$beta.exposure,
                                   bxse = dat_harm_t2dadjbmi_carb$se.exposure,
                                   by = dat_harm_t2dadjbmi_carb$beta.outcome,
                                   byse = dat_harm_t2dadjbmi_carb$se.outcome)
#Run MR
MRAll_t2dadjbmi_carb <- mr_allmethods(MRInput_t2dadjbmi_carb, method = "all")
res_t2dadjbmi_carb<- mr(dat_harm_t2dadjbmi_carb)
sink("YOUR_PATHWAY_DIRECTORY/MRAll_t2dadjbmi_carb.txt")
MRAll_t2dadjbmi_carb
generate_odds_ratios(res_t2dadjbmi_carb)
sink()
p_t2dadjbmi_carb<- mr_plot(mr_allmethods(MRInput_t2dadjbmi_carb, method = "main"))
p_t2dadjbmi_carb +
    xlab("Carbohydrates intake (E%)") +
    ylab("T2D adjusted for BMI (T2DadjBMI)")+
    ggsave(file = "YOUR_PATHWAY_DIRECTORY/MRAll_t2dadjbmi_carb_plot.pdf")

# Sensitivity analyses
het_comb_t2dadjbmi_carb<-mr_heterogeneity(dat_harm_t2dadjbmi_carb) #Heterogeneity Q and Q p-val
plt_comb_t2dadjbmi_carb<-mr_pleiotropy_test(dat_harm_t2dadjbmi_carb)  #Horizontal pleiotropy
sin_comb_t2dadjbmi_carb<-mr_singlesnp(dat_harm_t2dadjbmi_carb, all_method=c("mr_ivw", "mr_egger_regression"))
p_t2dadjbmi_carb <- mr_forest_plot(sin_comb_t2dadjbmi_carb)
p_t2dadjbmi_carb[[1]]
dir_t2dadjbmi_carb <- directionality_test(dat_harm_t2dadjbmi_carb)
leave_comb_t2dadjbmi_carb<- mr_leaveoneout(dat_harm_t2dadjbmi_carb, parameters = default_parameters())
p_leave_t2dadjbmi_carb <- mr_leaveoneout_plot(leave_comb_t2dadjbmi_carb)
p_leave_t2dadjbmi_carb[[1]]
#
mrpresso_t2dadjbmi_carb<- run_mr_presso(dat_harm_t2dadjbmi_carb, NbDistribution = 10000, SignifThreshold = 0.05)
saveRDS(mrpresso_t2dadjbmi_carb, file = "YOUR_PATHWAY_DIRECTORY/mr_presso/mrpresso_t2dadjbmi_carb.rds")

# T2D
# Harmonise
dat_harm_t2d_carb<- harmonise_data(
    exposure_dat = carb_exp_dat,
    outcome_dat = t2d_out_dat_carb)
#
mr_report(dat_harm_t2d_carb,output_path = "YOUR_PATHWAY_DIRECTORY/",
          output_type = "pdf",author = "Hugo P", study = "Two Sample MR",path = system.file("reports", package = "TwoSampleMR"))

# create input
MRInput_t2d_carb <- mr_input(bx = dat_harm_t2d_carb$beta.exposure,
                             bxse = dat_harm_t2d_carb$se.exposure,
                             by = dat_harm_t2d_carb$beta.outcome,
                             byse = dat_harm_t2d_carb$se.outcome)
#Run MR
MRAll_t2d_carb <- mr_allmethods(MRInput_t2d_carb, method = "all")
res_t2d_carb<- mr(dat_harm_t2d_carb)
sink("YOUR_PATHWAY_DIRECTORY/MRAll_t2d_carb.txt")
MRAll_t2d_carb
generate_odds_ratios(res_t2d_carb)
sink()
p_t2d_carb<- mr_plot(mr_allmethods(MRInput_t2d_carb, method = "main"))
p_t2d_carb +
    xlab("carbein intake (E%)") +
    ylab("T2D adjusted for BMI (T2D)") +
    ggsave(file = "YOUR_PATHWAY_DIRECTORY/MRAll_t2d_carb_plot.pdf")

# Sensitivity analyses
het_comb_t2d_carb<-mr_heterogeneity(dat_harm_t2d_carb) #Heterogeneity Q and Q p-val
plt_comb_t2d_carb<-mr_pleiotropy_test(dat_harm_t2d_carb)  #Horizontal pleiotropy
sin_comb_t2d_carb<-mr_singlesnp(dat_harm_t2d_carb, all_method=c("mr_ivw", "mr_egger_regression"))
p_t2d_carb <- mr_forest_plot(sin_comb_t2d_carb)
p_t2d_carb[[1]]
dir_t2d_carb <- directionality_test(dat_harm_t2d_carb)
leave_comb_t2d_carb<- mr_leaveoneout(dat_harm_t2d_carb, parameters = default_parameters())
p_leave_t2d_carb <- mr_leaveoneout_plot(leave_comb_t2d_carb)
p_leave_t2d_carb[[1]]
#
mrpresso_t2d_carb<- run_mr_presso(dat_harm_t2d_carb, NbDistribution = 10000, SignifThreshold = 0.05)
saveRDS(mrpresso_t2d_carb, file = "YOUR_PATHWAY_DIRECTORY/mrpresso_t2d_carb.rds")

# Fasting glucose
# Harmonise
dat_harm_fg_carb<- harmonise_data(
    exposure_dat = carb_exp_dat,
    outcome_dat = fg_out_dat_carb)
#
mr_report(dat_harm_fg_carb,output_path = "YOUR_PATHWAY_DIRECTORY/",
          output_type = "pdf",author = "Hugo P", study = "Two Sample MR",path = system.file("reports", package = "TwoSampleMR"))

# create input
MRInput_fg_carb <- mr_input(bx = dat_harm_fg_carb$beta.exposure,
                            bxse = dat_harm_fg_carb$se.exposure,
                            by = dat_harm_fg_carb$beta.outcome,
                            byse = dat_harm_fg_carb$se.outcome)
#Run MR
MRAll_fg_carb <- mr_allmethods(MRInput_fg_carb, method = "all")
sink("YOUR_PATHWAY_DIRECTORY/MRAll_fg_carb.txt")
MRAll_fg_carb
sink()
p_fg_carb<- mr_plot(mr_allmethods(MRInput_fg_carb, method = "main"))
p_fg_carb +
    xlab("Carbohydrates intake (E%)") +
    ylab("Fasting glucose (FG)")+
    ggsave(file = "YOUR_PATHWAY_DIRECTORY/MRAll_fg_carb_plot.pdf")

# Sensitivity analyses
het_comb_fg_carb<-mr_heterogeneity(dat_harm_fg_carb) #Heterogeneity Q and Q p-val
plt_comb_fg_carb<-mr_pleiotropy_test(dat_harm_fg_carb)  #Horizontal pleiotropy
sin_comb_fg_carb<-mr_singlesnp(dat_harm_fg_carb, all_method=c("mr_ivw", "mr_egger_regression"))
p_fg_carb <- mr_forest_plot(sin_comb_fg_carb)
p_fg_carb[[1]]
dir_fg_carb <- directionality_test(dat_harm_fg_carb)
leave_comb_fg_carb<- mr_leaveoneout(dat_harm_fg_carb, parameters = default_parameters())
p_leave_fg_carb <- mr_leaveoneout_plot(leave_comb_fg_carb)
p_leave_fg_carb[[1]]
#
mrpresso_fg_carb<- run_mr_presso(dat_harm_fg_carb, NbDistribution = 10000, SignifThreshold = 0.05)
saveRDS(mrpresso_fg_carb, file = "YOUR_PATHWAY_DIRECTORY/mrpresso_fg_carb.rds")
#  Not enough instrumental variables

# Two Hour glucose
# Harmonise
dat_harm_twohg_carb<- harmonise_data(
    exposure_dat = carb_exp_dat,
    outcome_dat = twohg_out_dat_carb)
#
mr_report(dat_harm_twohg_carb,output_path = "YOUR_PATHWAY_DIRECTORY/",
          output_type = "pdf",author = "Hugo P", study = "Two Sample MR",path = system.file("reports", package = "TwoSampleMR"))

# create input
MRInput_twohg_carb <- mr_input(bx = dat_harm_twohg_carb$beta.exposure,
                               bxse = dat_harm_twohg_carb$se.exposure,
                               by = dat_harm_twohg_carb$beta.outcome,
                               byse = dat_harm_twohg_carb$se.outcome)
#Run MR
MRAll_twohg_carb <- mr_allmethods(MRInput_twohg_carb, method = "all")
sink("YOUR_PATHWAY_DIRECTORY/MRAll_twohg_carb.txt")
MRAll_twohg_carb
sink()
p_twohg_carb<- mr_plot(mr_allmethods(MRInput_twohg_carb, method = "main"))
p_twohg_carb +
    xlab("Carbohydrates intake (E%)") +
    ylab("Two-hour glucose (2hG)")+
    ggsave(file = "YOUR_PATHWAY_DIRECTORY/MRAll_twohg_carb_plot.pdf")

# Sensitivity analyses
het_comb_twohg_carb<-mr_heterogeneity(dat_harm_twohg_carb) #Heterogeneity Q and Q p-val
plt_comb_twohg_carb<-mr_pleiotropy_test(dat_harm_twohg_carb)  #Horizontal pleiotropy
sin_comb_twohg_carb<-mr_singlesnp(dat_harm_twohg_carb, all_method=c("mr_ivw", "mr_egger_regression"))
p_twohg_carb <- mr_forest_plot(sin_comb_twohg_carb)
p_twohg_carb[[1]]
dir_twohg_carb <- directionality_test(dat_harm_twohg_carb)
leave_comb_twohg_carb<- mr_leaveoneout(dat_harm_twohg_carb, parameters = default_parameters())
p_leave_twohg_carb <- mr_leaveoneout_plot(leave_comb_twohg_carb)
p_leave_twohg_carb[[1]]
#
mrpresso_twohg_carb<- run_mr_presso(dat_harm_twohg_carb, NbDistribution = 10000, SignifThreshold = 0.05)
saveRDS(mrpresso_twohg_carb, file = "YOUR_PATHWAY_DIRECTORY/mrpresso_twohg_carb.rds")
#  Not enough instrumental variables

# Stroke
# Harmonise
dat_harm_stroke_carb<- harmonise_data(
    exposure_dat = carb_exp_dat,
    outcome_dat = stroke_out_dat_carb)
#
mr_report(dat_harm_stroke_carb,output_path = "YOUR_PATHWAY_DIRECTORY/",
          output_type = "pdf",author = "Hugo P", study = "Two Sample MR",path = system.file("reports", package = "TwoSampleMR"))

# create input
MRInput_stroke_carb <- mr_input(bx = dat_harm_stroke_carb$beta.exposure,
                                bxse = dat_harm_stroke_carb$se.exposure,
                                by = dat_harm_stroke_carb$beta.outcome,
                                byse = dat_harm_stroke_carb$se.outcome)
#Run MR
MRAll_stroke_carb <- mr_allmethods(MRInput_stroke_carb, method = "all") #Method requires data on >2 variants
res_stroke_carb<- mr(dat_harm_stroke_carb)
sink("YOUR_PATHWAY_DIRECTORY/MRAll_stroke_carb.txt")
MRAll_stroke_carb
generate_odds_ratios(res_stroke_carb)
sink()
p_stroke_carb<- mr_plot(mr_allmethods(MRInput_stroke_carb, method = "main"))
p_stroke_carb +
    xlab("Carbohydrates intake (E%)") +
    ylab("Stroke")+
    ggsave(file = "YOUR_PATHWAY_DIRECTORY/MRAll_stroke_carb_plot.pdf")

# Sensitivity analyses
het_comb_stroke_carb<-mr_heterogeneity(dat_harm_stroke_carb) #Heterogeneity Q and Q p-val
plt_comb_stroke_carb<-mr_pleiotropy_test(dat_harm_stroke_carb)  #Horizontal pleiotropy
sin_comb_stroke_carb<-mr_singlesnp(dat_harm_stroke_carb, all_method=c("mr_ivw", "mr_egger_regression"))
p_stroke_carb <- mr_forest_plot(sin_comb_stroke_carb)
p_stroke_carb[[1]]
dir_stroke_carb <- directionality_test(dat_harm_stroke_carb)
leave_comb_stroke_carb<- mr_leaveoneout(dat_harm_stroke_carb, parameters = default_parameters())
p_leave_stroke_carb <- mr_leaveoneout_plot(leave_comb_stroke_carb)
p_leave_stroke_carb[[1]]
#
mrpresso_stroke_carb<- run_mr_presso(dat_harm_stroke_carb, NbDistribution = 10000, SignifThreshold = 0.05)
saveRDS(mrpresso_stroke_carb, file = "YOUR_PATHWAY_DIRECTORY/mrpresso_stroke_carb.rds")
#  Not enough instrumental variables

# HDL
# Harmonise
dat_harm_hdl_carb<- harmonise_data(
    exposure_dat = carb_exp_dat,
    outcome_dat = hdl_out_dat_carb)
#
mr_report(dat_harm_hdl_carb,output_path = "YOUR_PATHWAY_DIRECTORY/",
          output_type = "pdf",author = "Hugo P", study = "Two Sample MR",path = system.file("reports", package = "TwoSampleMR"))

# create input
MRInput_hdl_carb <- mr_input(bx = dat_harm_hdl_carb$beta.exposure,
                             bxse = dat_harm_hdl_carb$se.exposure,
                             by = dat_harm_hdl_carb$beta.outcome,
                             byse = dat_harm_hdl_carb$se.outcome)
#Run MR
MRAll_hdl_carb <- mr_allmethods(MRInput_hdl_carb, method = "all")
sink("YOUR_PATHWAY_DIRECTORY/MRAll_hdl_carb.txt")
MRAll_hdl_carb
sink()
p_hdl_carb<- mr_plot(mr_allmethods(MRInput_hdl_carb, method = "main"))
p_hdl_carb +
    xlab("carbein intake (E%)") +
    ylab("hdl")+
    ggsave(file = "YOUR_PATHWAY_DIRECTORY/MRAll_hdl_carb_plot.pdf")

# Sensitivity analyses
het_comb_hdl_carb<-mr_heterogeneity(dat_harm_hdl_carb) #Heterogeneity Q and Q p-val
plt_comb_hdl_carb<-mr_pleiotropy_test(dat_harm_hdl_carb)  #Horizontal pleiotropy
sin_comb_hdl_carb<-mr_singlesnp(dat_harm_hdl_carb, all_method=c("mr_ivw", "mr_egger_regression"))
p_hdl_carb <- mr_forest_plot(sin_comb_hdl_carb)
p_hdl_carb[[1]]
dir_hdl_carb <- directionality_test(dat_harm_hdl_carb)
leave_comb_hdl_carb<- mr_leaveoneout(dat_harm_hdl_carb, parameters = default_parameters())
p_leave_hdl_carb <- mr_leaveoneout_plot(leave_comb_hdl_carb)
p_leave_hdl_carb[[1]]
#
mrpresso_hdl_carb<- run_mr_presso(dat_harm_hdl_carb, NbDistribution = 10000, SignifThreshold = 0.05)
saveRDS(mrpresso_hdl_carb, file = "YOUR_PATHWAY_DIRECTORY/mrpresso_hdl_carb.rds")

# LDL
# Harmonise
dat_harm_ldl_carb<- harmonise_data(
    exposure_dat = carb_exp_dat,
    outcome_dat = ldl_out_dat_carb)
#
mr_report(dat_harm_ldl_carb,output_path = "YOUR_PATHWAY_DIRECTORY/",
          output_type = "pdf",author = "Hugo P", study = "Two Sample MR",path = system.file("reports", package = "TwoSampleMR"))

# create input
MRInput_ldl_carb <- mr_input(bx = dat_harm_ldl_carb$beta.exposure,
                             bxse = dat_harm_ldl_carb$se.exposure,
                             by = dat_harm_ldl_carb$beta.outcome,
                             byse = dat_harm_ldl_carb$se.outcome)
#Run MR
MRAll_ldl_carb <- mr_allmethods(MRInput_ldl_carb, method = "all")
sink("YOUR_PATHWAY_DIRECTORY/MRAll_ldl_carb.txt")
MRAll_ldl_carb
sink()
p_ldl_carb<- mr_plot(mr_allmethods(MRInput_ldl_carb, method = "main"))
p_ldl_carb +
    xlab("carbein intake (E%)") +
    ylab("ldl")+
    ggsave(file = "YOUR_PATHWAY_DIRECTORY/MRAll_ldl_carb_plot.pdf")

# Sensitivity analyses
het_comb_ldl_carb<-mr_heterogeneity(dat_harm_ldl_carb) #Heterogeneity Q and Q p-val
plt_comb_ldl_carb<-mr_pleiotropy_test(dat_harm_ldl_carb)  #Horizontal pleiotropy
sin_comb_ldl_carb<-mr_singlesnp(dat_harm_ldl_carb, all_method=c("mr_ivw", "mr_egger_regression"))
p_ldl_carb <- mr_forest_plot(sin_comb_ldl_carb)
p_ldl_carb[[1]]
dir_ldl_carb <- directionality_test(dat_harm_ldl_carb)
leave_comb_ldl_carb<- mr_leaveoneout(dat_harm_ldl_carb, parameters = default_parameters())
p_leave_ldl_carb <- mr_leaveoneout_plot(leave_comb_ldl_carb)
p_leave_ldl_carb[[1]]
#
mrpresso_ldl_carb<- run_mr_presso(dat_harm_ldl_carb, NbDistribution = 10000, SignifThreshold = 0.05)
saveRDS(mrpresso_ldl_carb, file = "YOUR_PATHWAY_DIRECTORY/mrpresso_ldl_carb.rds")

# CHOLESTEROL
# Harmonise
dat_harm_skol_carb<- harmonise_data(
    exposure_dat = carb_exp_dat,
    outcome_dat = skol_out_dat_carb)
#
mr_report(dat_harm_skol_carb,output_path = "YOUR_PATHWAY_DIRECTORY/",
          output_type = "pdf",author = "Hugo P", study = "Two Sample MR",path = system.file("reports", package = "TwoSampleMR"))

# create input
MRInput_skol_carb <- mr_input(bx = dat_harm_skol_carb$beta.exposure,
                              bxse = dat_harm_skol_carb$se.exposure,
                              by = dat_harm_skol_carb$beta.outcome,
                              byse = dat_harm_skol_carb$se.outcome)
#Run MR
MRAll_skol_carb <- mr_allmethods(MRInput_skol_carb, method = "all")
sink("YOUR_PATHWAY_DIRECTORY/MRAll_skol_carb.txt")
MRAll_skol_carb
sink()
p_skol_carb<- mr_plot(mr_allmethods(MRInput_skol_carb, method = "main"))
p_skol_carb +
    xlab("carbein intake (E%)") +
    ylab("skol")+
    ggsave(file = "YOUR_PATHWAY_DIRECTORY/MRAll_skol_carb_plot.pdf")

# Sensitivity analyses
het_comb_skol_carb<-mr_heterogeneity(dat_harm_skol_carb) #Heterogeneity Q and Q p-val
plt_comb_skol_carb<-mr_pleiotropy_test(dat_harm_skol_carb)  #Horizontal pleiotropy
sin_comb_skol_carb<-mr_singlesnp(dat_harm_skol_carb, all_method=c("mr_ivw", "mr_egger_regression"))
p_skol_carb <- mr_forest_plot(sin_comb_skol_carb)
p_skol_carb[[1]]
dir_skol_carb <- directionality_test(dat_harm_skol_carb)
leave_comb_skol_carb<- mr_leaveoneout(dat_harm_skol_carb, parameters = default_parameters())
p_leave_skol_carb <- mr_leaveoneout_plot(leave_comb_skol_carb)
p_leave_skol_carb[[1]]
#
mrpresso_skol_carb<- run_mr_presso(dat_harm_skol_carb, NbDistribution = 10000, SignifThreshold = 0.05)
saveRDS(mrpresso_skol_carb, file = "YOUR_PATHWAY_DIRECTORY/mrpresso_skol_carb.rds")

# TRIGLYCERIDES
# Harmonise
dat_harm_stg_carb<- harmonise_data(
    exposure_dat = carb_exp_dat,
    outcome_dat = stg_out_dat_carb)
#
mr_report(dat_harm_stg_carb,output_path = "YOUR_PATHWAY_DIRECTORY/",
          output_type = "pdf",author = "Hugo P", study = "Two Sample MR",path = system.file("reports", package = "TwoSampleMR"))

# create input
MRInput_stg_carb <- mr_input(bx = dat_harm_stg_carb$beta.exposure,
                             bxse = dat_harm_stg_carb$se.exposure,
                             by = dat_harm_stg_carb$beta.outcome,
                             byse = dat_harm_stg_carb$se.outcome)
#Run MR
MRAll_stg_carb <- mr_allmethods(MRInput_stg_carb, method = "all")
sink("YOUR_PATHWAY_DIRECTORY/MRAll_stg_carb.txt")
MRAll_stg_carb
sink()
p_stg_carb<- mr_plot(mr_allmethods(MRInput_stg_carb, method = "main"))
p_stg_carb +
    xlab("carbein intake (E%)") +
    ylab("stg")+
    ggsave(file = "YOUR_PATHWAY_DIRECTORY/MRAll_stg_carb_plot.pdf")

# Sensitivity analyses
het_comb_stg_carb<-mr_heterogeneity(dat_harm_stg_carb) #Heterogeneity Q and Q p-val
plt_comb_stg_carb<-mr_pleiotropy_test(dat_harm_stg_carb)  #Horizontal pleiotropy
sin_comb_stg_carb<-mr_singlesnp(dat_harm_stg_carb, all_method=c("mr_ivw", "mr_egger_regression"))
p_stg_carb <- mr_forest_plot(sin_comb_stg_carb)
p_stg_carb[[1]]
dir_stg_carb <- directionality_test(dat_harm_stg_carb)
leave_comb_stg_carb<- mr_leaveoneout(dat_harm_stg_carb, parameters = default_parameters())
p_leave_stg_carb <- mr_leaveoneout_plot(leave_comb_stg_carb)
p_leave_stg_carb[[1]]
#
mrpresso_stg_carb<- run_mr_presso(dat_harm_stg_carb, NbDistribution = 10000, SignifThreshold = 0.05)
saveRDS(mrpresso_stg_carb, file = "YOUR_PATHWAY_DIRECTORY/mrpresso_stg_carb.rds")

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#Read outcomes for sugars
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
chd_out_dat_sugar <- extract_outcome_data(
    snps = sugar_exp_dat$SNP,
    outcomes = 'ieu-a-7'
)

t2dadjbmi_out_dat_sugar <- extract_outcome_data(
    snps = sugar_exp_dat$SNP,
    outcomes = 'ebi-a-GCST007518' #GCST007515
)
t2d_out_dat_sugar <- extract_outcome_data(
    snps = sugar_exp_dat$SNP,
    outcomes = 'ebi-a-GCST007517'
)
fg_out_dat_sugar <- extract_outcome_data(
    snps = sugar_exp_dat$SNP,
    outcomes = 'ieu-b-113'
)

twohg_out_dat_sugar <- extract_outcome_data(
    snps = sugar_exp_dat$SNP,
    outcomes = 'ebi-a-GCST000569'
)
stroke_out_dat_sugar <- extract_outcome_data(
    snps = prot_exp_dat$SNP,
    outcomes = 'ebi-a-GCST006906'
)
hdl_out_dat_sugar <- extract_outcome_data(
    snps = sugar_exp_dat$SNP,
    outcomes = 'ieu-b-109'
)
ldl_out_dat_sugar <- extract_outcome_data(
    snps = sugar_exp_dat$SNP,
    outcomes = 'ieu-b-110'
)
skol_out_dat_sugar <- extract_outcome_data(
    snps = sugar_exp_dat$SNP,
    outcomes = 'met-d-Total_C'
)
stg_out_dat_sugar <- extract_outcome_data(
    snps = sugar_exp_dat$SNP,
    outcomes = 'ieu-b-111'
)

# CHD
# Harmonise
dat_harm_chd_sugar<- harmonise_data(
    exposure_dat = sugar_exp_dat,
    outcome_dat = chd_out_dat_sugar)
#
mr_report(dat_harm_chd_sugar,output_path = "YOUR_PATHWAY_DIRECTORY/",
          output_type = "pdf",author = "Hugo P", study = "Two Sample MR",path = system.file("reports", package = "TwoSampleMR"))

# create input
MRInput_chd_sugar <- mr_input(bx = dat_harm_chd_sugar$beta.exposure,
                              bxse = dat_harm_chd_sugar$se.exposure,
                              by = dat_harm_chd_sugar$beta.outcome,
                              byse = dat_harm_chd_sugar$se.outcome)
#Run MR
MRAll_chd_sugar <- mr_allmethods(MRInput_chd_sugar, method = "all")
res_chd_sugar<- mr(dat_harm_chd_sugar)
sink("YOUR_PATHWAY_DIRECTORY/MRAll_chd_sugar.txt")
MRAll_chd_sugar
generate_odds_ratios(res_chd_sugar)
sink()
p_chd_sugar<- mr_plot(mr_allmethods(MRInput_chd_sugar, method = "main"))
p_chd_sugar +
    xlab("Sugar intake (E%)") +
    ylab("Coronary heart disease (CHD)")+
    ggsave(file = "YOUR_PATHWAY_DIRECTORY/MRAll_chd_sugar_plot.pdf")

# Sensitivity analyses
het_comb_chd_sugar<-mr_heterogeneity(dat_harm_chd_sugar) #Heterogeneity Q and Q p-val
plt_comb_chd_sugar<-mr_pleiotropy_test(dat_harm_chd_sugar)  #Horizontal pleiotropy
sin_comb_chd_sugar<-mr_singlesnp(dat_harm_chd_sugar, all_method=c("mr_ivw", "mr_egger_regression"))
p_chd_sugar <- mr_forest_plot(sin_comb_chd_sugar)
p_chd_sugar[[1]]
dir_chd_sugar <- directionality_test(dat_harm_chd_sugar)
leave_comb_chd_sugar<- mr_leaveoneout(dat_harm_chd_sugar, parameters = default_parameters())
p_leave_chd_sugar <- mr_leaveoneout_plot(leave_comb_chd_sugar)
p_leave_chd_sugar[[1]]
#
mrpresso_chd_sugar<- run_mr_presso(dat_harm_chd_sugar, NbDistribution = 10000, SignifThreshold = 0.05)
saveRDS(mrpresso_chd_sugar, file = "YOUR_PATHWAY_DIRECTORY/mrpresso_chd_sugar.rds")

# t2dadjbmi
# Harmonise
dat_harm_t2dadjbmi_sugar<- harmonise_data(
    exposure_dat = sugar_exp_dat,
    outcome_dat = t2dadjbmi_out_dat_sugar)
#
mr_report(dat_harm_t2dadjbmi_sugar,output_path = "YOUR_PATHWAY_DIRECTORY/",
          output_type = "pdf",author = "Hugo P", study = "Two Sample MR",path = system.file("reports", package = "TwoSampleMR"))

# create input
MRInput_t2dadjbmi_sugar <- mr_input(bx = dat_harm_t2dadjbmi_sugar$beta.exposure,
                                    bxse = dat_harm_t2dadjbmi_sugar$se.exposure,
                                    by = dat_harm_t2dadjbmi_sugar$beta.outcome,
                                    byse = dat_harm_t2dadjbmi_sugar$se.outcome)
#Run MR
MRAll_t2dadjbmi_sugar <- mr_allmethods(MRInput_t2dadjbmi_sugar, method = "all")
res_t2dadjbmi_sugar<- mr(dat_harm_t2dadjbmi_sugar)
sink("YOUR_PATHWAY_DIRECTORY/MRAll_t2dadjbmi_sugar.txt")
MRAll_t2dadjbmi_sugar
generate_odds_ratios(res_t2dadjbmi_sugar)
sink()
#Method requires data on >2 variants
p_t2dadjbmi_sugar<- mr_plot(mr_allmethods(MRInput_t2dadjbmi_sugar, method = "main"))
p_t2dadjbmi_sugar +
    xlab("Sugar intake (E%)") +
    ylab("T2D adjusted for BMI (T2DadjBMI)")+
    ggsave(file = "YOUR_PATHWAY_DIRECTORY/MRAll_t2dadjbmi_sugar_plot.pdf")

# Sensitivity analyses
het_comb_t2dadjbmi_sugar<-mr_heterogeneity(dat_harm_t2dadjbmi_sugar) #Heterogeneity Q and Q p-val
plt_comb_t2dadjbmi_sugar<-mr_pleiotropy_test(dat_harm_t2dadjbmi_sugar)  #Horizontal pleiotropy
sin_comb_t2dadjbmi_sugar<-mr_singlesnp(dat_harm_t2dadjbmi_sugar, all_method=c("mr_ivw", "mr_egger_regression"))
p_t2dadjbmi_sugar <- mr_forest_plot(sin_comb_t2dadjbmi_sugar)
p_t2dadjbmi_sugar[[1]]
dir_t2dadjbmi_sugar <- directionality_test(dat_harm_t2dadjbmi_sugar)
leave_comb_t2dadjbmi_sugar<- mr_leaveoneout(dat_harm_t2dadjbmi_sugar, parameters = default_parameters())
p_leave_t2dadjbmi_sugar <- mr_leaveoneout_plot(leave_comb_t2dadjbmi_sugar)
p_leave_t2dadjbmi_sugar[[1]]
#
mrpresso_t2dadjbmi_sugar<- run_mr_presso(dat_harm_t2dadjbmi_sugar, NbDistribution = 10000, SignifThreshold = 0.05)
#Not enough SNPs available
saveRDS(mrpresso_t2dadjbmi_sugar, file = "YOUR_PATHWAY_DIRECTORY/mrpresso_t2dadjbmi_sugar.rds")

# t2d
# Harmonise
dat_harm_t2d_sugar<- harmonise_data(
    exposure_dat = sugar_exp_dat,
    outcome_dat = t2d_out_dat_sugar)
#
mr_report(dat_harm_t2d_sugar,output_path = "YOUR_PATHWAY_DIRECTORY/",
          output_type = "pdf",author = "Hugo P", study = "Two Sample MR",path = system.file("reports", package = "TwoSampleMR"))

# create input
MRInput_t2d_sugar <- mr_input(bx = dat_harm_t2d_sugar$beta.exposure,
                              bxse = dat_harm_t2d_sugar$se.exposure,
                              by = dat_harm_t2d_sugar$beta.outcome,
                              byse = dat_harm_t2d_sugar$se.outcome)
#Run MR
MRAll_t2d_sugar <- mr_allmethods(MRInput_t2d_sugar, method = "all")
res_t2d_sugar<- mr(dat_harm_t2d_sugar)
sink("YOUR_PATHWAY_DIRECTORY/MRAll_t2d_sugar.txt")
MRAll_t2d_sugar
generate_odds_ratios(res_t2d_sugar)
sink()
p_t2d_sugar<- mr_plot(mr_allmethods(MRInput_t2d_sugar, method = "main"))
p_t2d_sugar +
    xlab("sugarein intake (E%)") +
    ylab("T2D adjusted for BMI (T2D)") +
    ggsave(file = "YOUR_PATHWAY_DIRECTORY/MRAll_t2d_sugar_plot.pdf")

# Sensitivity analyses
het_comb_t2d_sugar<-mr_heterogeneity(dat_harm_t2d_sugar) #Heterogeneity Q and Q p-val
plt_comb_t2d_sugar<-mr_pleiotropy_test(dat_harm_t2d_sugar)  #Horizontal pleiotropy
sin_comb_t2d_sugar<-mr_singlesnp(dat_harm_t2d_sugar, all_method=c("mr_ivw", "mr_egger_regression"))
p_t2d_sugar <- mr_forest_plot(sin_comb_t2d_sugar)
p_t2d_sugar[[1]]
dir_t2d_sugar <- directionality_test(dat_harm_t2d_sugar)
leave_comb_t2d_sugar<- mr_leaveoneout(dat_harm_t2d_sugar, parameters = default_parameters())
p_leave_t2d_sugar <- mr_leaveoneout_plot(leave_comb_t2d_sugar)
p_leave_t2d_sugar[[1]]
#
mrpresso_t2d_sugar<- run_mr_presso(dat_harm_t2d_sugar, NbDistribution = 10000, SignifThreshold = 0.05)
saveRDS(mrpresso_t2d_sugar, file = "YOUR_PATHWAY_DIRECTORY/mrpresso_t2d_sugar.rds")
# Not enough intrumental variables

# Fasting glucose
# Harmonise
dat_harm_fg_sugar<- harmonise_data(
    exposure_dat = sugar_exp_dat,
    outcome_dat = fg_out_dat_sugar)
#
mr_report(dat_harm_fg_sugar,output_path = "YOUR_PATHWAY_DIRECTORY/",
          output_type = "pdf",author = "Hugo P", study = "Two Sample MR",path = system.file("reports", package = "TwoSampleMR"))

# create input
MRInput_fg_sugar <- mr_input(bx = dat_harm_fg_sugar$beta.exposure,
                             bxse = dat_harm_fg_sugar$se.exposure,
                             by = dat_harm_fg_sugar$beta.outcome,
                             byse = dat_harm_fg_sugar$se.outcome)
#Run MR
MRAll_fg_sugar <- mr_allmethods(MRInput_fg_sugar, method = "all")
sink("YOUR_PATHWAY_DIRECTORY/MRAll_fg_sugar.txt")
MRAll_fg_sugar
sink()
p_fg_sugar<- mr_plot(mr_allmethods(MRInput_fg_sugar, method = "main"))
p_fg_sugar +
    xlab("Sugar intake (E%)") +
    ylab("T2D adjusted for BMI (fg)")+
    ggsave(file = "YOUR_PATHWAY_DIRECTORY/MRAll_fg_sugar_plot.pdf")

# Sensitivity analyses
het_comb_fg_sugar<-mr_heterogeneity(dat_harm_fg_sugar) #Heterogeneity Q and Q p-val
plt_comb_fg_sugar<-mr_pleiotropy_test(dat_harm_fg_sugar)  #Horizontal pleiotropy
sin_comb_fg_sugar<-mr_singlesnp(dat_harm_fg_sugar, all_method=c("mr_ivw", "mr_egger_regression"))
p_fg_sugar <- mr_forest_plot(sin_comb_fg_sugar)
p_fg_sugar[[1]]
dir_fg_sugar <- directionality_test(dat_harm_fg_sugar)
leave_comb_fg_sugar<- mr_leaveoneout(dat_harm_fg_sugar, parameters = default_parameters())
p_leave_fg_sugar <- mr_leaveoneout_plot(leave_comb_fg_sugar)
p_leave_fg_sugar[[1]]
#
mrpresso_fg_sugar<- run_mr_presso(dat_harm_fg_sugar, NbDistribution = 10000, SignifThreshold = 0.05)
saveRDS(mrpresso_fg_sugar, file = "mrpresso_fg_sugar.rds")
#  Not enough instrumental variables

# Two Hour glucose
# Harmonise
dat_harm_twohg_sugar<- harmonise_data(
    exposure_dat = sugar_exp_dat,
    outcome_dat = twohg_out_dat_sugar)
#
mr_report(dat_harm_twohg_sugar,output_path = "YOUR_PATHWAY_DIRECTORY/",
          output_type = "pdf",author = "Hugo P", study = "Two Sample MR",path = system.file("reports", package = "TwoSampleMR"))

# create input
MRInput_twohg_sugar <- mr_input(bx = dat_harm_twohg_sugar$beta.exposure,
                                bxse = dat_harm_twohg_sugar$se.exposure,
                                by = dat_harm_twohg_sugar$beta.outcome,
                                byse = dat_harm_twohg_sugar$se.outcome)
#Run MR
MRAll_twohg_sugar <- mr_allmethods(MRInput_twohg_sugar, method = "all")
sink("YOUR_PATHWAY_DIRECTORY/MRAll_twohg_sugar.txt")
MRAll_twohg_sugar
sink()
p_twohg_sugar<- mr_plot(mr_allmethods(MRInput_twohg_sugar, method = "main"))
p_twohg_sugar +
    xlab("Sugar intake (E%)") +
    ylab("Two-hour glucose (2hG)")+
    ggsave(file = "YOUR_PATHWAY_DIRECTORY/MRAll_twohg_sugar_plot.pdf")

# Sensitivity analyses
het_comb_twohg_sugar<-mr_heterogeneity(dat_harm_twohg_sugar) #Heterogeneity Q and Q p-val
plt_comb_twohg_sugar<-mr_pleiotropy_test(dat_harm_twohg_sugar)  #Horizontal pleiotropy
sin_comb_twohg_sugar<-mr_singlesnp(dat_harm_twohg_sugar, all_method=c("mr_ivw", "mr_egger_regression"))
p_twohg_sugar <- mr_forest_plot(sin_comb_twohg_sugar)
p_twohg_sugar[[1]]
dir_twohg_sugar <- directionality_test(dat_harm_twohg_sugar)
leave_comb_twohg_sugar<- mr_leaveoneout(dat_harm_twohg_sugar, parameters = default_parameters())
p_leave_twohg_sugar <- mr_leaveoneout_plot(leave_comb_twohg_sugar)
p_leave_twohg_sugar[[1]]
#
mrpresso_twohg_sugar<- run_mr_presso(dat_harm_twohg_sugar, NbDistribution = 10000, SignifThreshold = 0.05)
saveRDS(mrpresso_twohg_sugar, file = "mrpresso_twohg_sugar.rds")
#  Not enough instrumental variables

# Stroke
# Harmonise
dat_harm_stroke_sugar<- harmonise_data(
    exposure_dat = sugar_exp_dat,
    outcome_dat = stroke_out_dat_sugar)
#
mr_report(dat_harm_stroke_sugar,output_path = "YOUR_PATHWAY_DIRECTORY/",
          output_type = "pdf",author = "Hugo P", study = "Two Sample MR",path = system.file("reports", package = "TwoSampleMR"))

# create input
MRInput_stroke_sugar <- mr_input(bx = dat_harm_stroke_sugar$beta.exposure,
                                 bxse = dat_harm_stroke_sugar$se.exposure,
                                 by = dat_harm_stroke_sugar$beta.outcome,
                                 byse = dat_harm_stroke_sugar$se.outcome)
#Run MR
MRAll_stroke_sugar <- mr_allmethods(MRInput_stroke_sugar, method = "all") #Method requires data on >2 variants
res_stroke_sugar<- mr(dat_harm_stroke_sugar)
sink("YOUR_PATHWAY_DIRECTORY/MRAll_stroke_sugar.txt")
MRAll_stroke_sugar
generate_odds_ratios(res_stroke_sugar)
sink()
p_stroke_sugar<- mr_plot(mr_allmethods(MRInput_stroke_sugar, method = "main"))
p_stroke_sugar +
    xlab("Sugar intake (E%)") +
    ylab("Stroke")+
    ggsave(file = "YOUR_PATHWAY_DIRECTORY/MRAll_stroke_sugar_plot.pdf")

# Sensitivity analyses
het_comb_stroke_sugar<-mr_heterogeneity(dat_harm_stroke_sugar) #Heterogeneity Q and Q p-val
plt_comb_stroke_sugar<-mr_pleiotropy_test(dat_harm_stroke_sugar)  #Horizontal pleiotropy
sin_comb_stroke_sugar<-mr_singlesnp(dat_harm_stroke_sugar, all_method=c("mr_ivw", "mr_egger_regression"))
p_stroke_sugar <- mr_forest_plot(sin_comb_stroke_sugar)
p_stroke_sugar[[1]]
dir_stroke_sugar <- directionality_test(dat_harm_stroke_sugar)
leave_comb_stroke_sugar<- mr_leaveoneout(dat_harm_stroke_sugar, parameters = default_parameters())
p_leave_stroke_sugar <- mr_leaveoneout_plot(leave_comb_stroke_sugar)
p_leave_stroke_sugar[[1]]
#
mrpresso_stroke_sugar<- run_mr_presso(dat_harm_stroke_sugar, NbDistribution = 10000, SignifThreshold = 0.05)
saveRDS(mrpresso_stroke_sugar, file = "YOUR_PATHWAY_DIRECTORY/mrpresso_stroke_sugar.rds")
#  Not enough instrumental variables

# HDL
# Harmonise
dat_harm_hdl_sugar<- harmonise_data(
    exposure_dat = sugar_exp_dat,
    outcome_dat = hdl_out_dat_sugar)
#
mr_report(dat_harm_hdl_sugar,output_path = "YOUR_PATHWAY_DIRECTORY/",
          output_type = "pdf",author = "Hugo P", study = "Two Sample MR",path = system.file("reports", package = "TwoSampleMR"))

# create input
MRInput_hdl_sugar <- mr_input(bx = dat_harm_hdl_sugar$beta.exposure,
                              bxse = dat_harm_hdl_sugar$se.exposure,
                              by = dat_harm_hdl_sugar$beta.outcome,
                              byse = dat_harm_hdl_sugar$se.outcome)
#Run MR
MRAll_hdl_sugar <- mr_allmethods(MRInput_hdl_sugar, method = "all")
sink("YOUR_PATHWAY_DIRECTORY/MRAll_hdl_sugar.txt")
MRAll_hdl_sugar
sink()
p_hdl_sugar<- mr_plot(mr_allmethods(MRInput_hdl_sugar, method = "main"))
p_hdl_sugar +
    xlab("sugarein intake (E%)") +
    ylab("hdl")+
    ggsave(file = "YOUR_PATHWAY_DIRECTORY/MRAll_hdl_sugar_plot.pdf")

# Sensitivity analyses
het_comb_hdl_sugar<-mr_heterogeneity(dat_harm_hdl_sugar) #Heterogeneity Q and Q p-val
plt_comb_hdl_sugar<-mr_pleiotropy_test(dat_harm_hdl_sugar)  #Horizontal pleiotropy
sin_comb_hdl_sugar<-mr_singlesnp(dat_harm_hdl_sugar, all_method=c("mr_ivw", "mr_egger_regression"))
p_hdl_sugar <- mr_forest_plot(sin_comb_hdl_sugar)
p_hdl_sugar[[1]]
dir_hdl_sugar <- directionality_test(dat_harm_hdl_sugar)
leave_comb_hdl_sugar<- mr_leaveoneout(dat_harm_hdl_sugar, parameters = default_parameters())
p_leave_hdl_sugar <- mr_leaveoneout_plot(leave_comb_hdl_sugar)
p_leave_hdl_sugar[[1]]
#
mrpresso_hdl_sugar<- run_mr_presso(dat_harm_hdl_sugar, NbDistribution = 10000, SignifThreshold = 0.05)
saveRDS(mrpresso_hdl_sugar, file = "YOUR_PATHWAY_DIRECTORY/mrpresso_hdl_sugar.rds")

# LDL
# Harmonise
dat_harm_ldl_sugar<- harmonise_data(
    exposure_dat = sugar_exp_dat,
    outcome_dat = ldl_out_dat_sugar)
#
mr_report(dat_harm_ldl_sugar,output_path = "YOUR_PATHWAY_DIRECTORY/",
          output_type = "pdf",author = "Hugo P", study = "Two Sample MR",path = system.file("reports", package = "TwoSampleMR"))

# create input
MRInput_ldl_sugar <- mr_input(bx = dat_harm_ldl_sugar$beta.exposure,
                              bxse = dat_harm_ldl_sugar$se.exposure,
                              by = dat_harm_ldl_sugar$beta.outcome,
                              byse = dat_harm_ldl_sugar$se.outcome)
#Run MR
MRAll_ldl_sugar <- mr_allmethods(MRInput_ldl_sugar, method = "all")
sink("YOUR_PATHWAY_DIRECTORY/MRAll_ldl_sugar.txt")
MRAll_ldl_sugar
sink()
p_ldl_sugar<- mr_plot(mr_allmethods(MRInput_ldl_sugar, method = "main"))
p_ldl_sugar +
    xlab("sugarein intake (E%)") +
    ylab("ldl")+
    ggsave(file = "YOUR_PATHWAY_DIRECTORY/MRAll_ldl_sugar_plot.pdf")

# Sensitivity analyses
het_comb_ldl_sugar<-mr_heterogeneity(dat_harm_ldl_sugar) #Heterogeneity Q and Q p-val
plt_comb_ldl_sugar<-mr_pleiotropy_test(dat_harm_ldl_sugar)  #Horizontal pleiotropy
sin_comb_ldl_sugar<-mr_singlesnp(dat_harm_ldl_sugar, all_method=c("mr_ivw", "mr_egger_regression"))
p_ldl_sugar <- mr_forest_plot(sin_comb_ldl_sugar)
p_ldl_sugar[[1]]
dir_ldl_sugar <- directionality_test(dat_harm_ldl_sugar)
leave_comb_ldl_sugar<- mr_leaveoneout(dat_harm_ldl_sugar, parameters = default_parameters())
p_leave_ldl_sugar <- mr_leaveoneout_plot(leave_comb_ldl_sugar)
p_leave_ldl_sugar[[1]]
#
mrpresso_ldl_sugar<- run_mr_presso(dat_harm_ldl_sugar, NbDistribution = 10000, SignifThreshold = 0.05)
saveRDS(mrpresso_ldl_sugar, file = "YOUR_PATHWAY_DIRECTORY/mrpresso_ldl_sugar.rds")

# CHOLESTEROL
# Harmonise
dat_harm_skol_sugar<- harmonise_data(
    exposure_dat = sugar_exp_dat,
    outcome_dat = skol_out_dat_sugar)
#
mr_report(dat_harm_skol_sugar,output_path = "YOUR_PATHWAY_DIRECTORY/",
          output_type = "pdf",author = "Hugo P", study = "Two Sample MR",path = system.file("reports", package = "TwoSampleMR"))

# create input
MRInput_skol_sugar <- mr_input(bx = dat_harm_skol_sugar$beta.exposure,
                               bxse = dat_harm_skol_sugar$se.exposure,
                               by = dat_harm_skol_sugar$beta.outcome,
                               byse = dat_harm_skol_sugar$se.outcome)
#Run MR
MRAll_skol_sugar <- mr_allmethods(MRInput_skol_sugar, method = "all")
sink("YOUR_PATHWAY_DIRECTORY/MRAll_skol_sugar.txt")
MRAll_skol_sugar
sink()
p_skol_sugar<- mr_plot(mr_allmethods(MRInput_skol_sugar, method = "main"))
p_skol_sugar +
    xlab("sugarein intake (E%)") +
    ylab("skol")+
    ggsave(file = "YOUR_PATHWAY_DIRECTORY/MRAll_skol_sugar_plot.pdf")

# Sensitivity analyses
het_comb_skol_sugar<-mr_heterogeneity(dat_harm_skol_sugar) #Heterogeneity Q and Q p-val
plt_comb_skol_sugar<-mr_pleiotropy_test(dat_harm_skol_sugar)  #Horizontal pleiotropy
sin_comb_skol_sugar<-mr_singlesnp(dat_harm_skol_sugar, all_method=c("mr_ivw", "mr_egger_regression"))
p_skol_sugar <- mr_forest_plot(sin_comb_skol_sugar)
p_skol_sugar[[1]]
dir_skol_sugar <- directionality_test(dat_harm_skol_sugar)
leave_comb_skol_sugar<- mr_leaveoneout(dat_harm_skol_sugar, parameters = default_parameters())
p_leave_skol_sugar <- mr_leaveoneout_plot(leave_comb_skol_sugar)
p_leave_skol_sugar[[1]]
#
mrpresso_skol_sugar<- run_mr_presso(dat_harm_skol_sugar, NbDistribution = 10000, SignifThreshold = 0.05)
saveRDS(mrpresso_skol_sugar, file = "YOUR_PATHWAY_DIRECTORY/mrpresso_skol_sugar.rds")

# TRIGLYCERIDES
# Harmonise
dat_harm_stg_sugar<- harmonise_data(
    exposure_dat = sugar_exp_dat,
    outcome_dat = stg_out_dat_sugar)
#
mr_report(dat_harm_stg_sugar,output_path = "YOUR_PATHWAY_DIRECTORY/",
          output_type = "pdf",author = "Hugo P", study = "Two Sample MR",path = system.file("reports", package = "TwoSampleMR"))

# create input
MRInput_stg_sugar <- mr_input(bx = dat_harm_stg_sugar$beta.exposure,
                              bxse = dat_harm_stg_sugar$se.exposure,
                              by = dat_harm_stg_sugar$beta.outcome,
                              byse = dat_harm_stg_sugar$se.outcome)
#Run MR
MRAll_stg_sugar <- mr_allmethods(MRInput_stg_sugar, method = "all")
sink("YOUR_PATHWAY_DIRECTORY/MRAll_stg_sugar.txt")
MRAll_stg_sugar
sink()
p_stg_sugar<- mr_plot(mr_allmethods(MRInput_stg_sugar, method = "main"))
p_stg_sugar +
    xlab("sugarein intake (E%)") +
    ylab("stg")+
    ggsave(file = "YOUR_PATHWAY_DIRECTORY/MRAll_stg_sugar_plot.pdf")

# Sensitivity analyses
het_comb_stg_sugar<-mr_heterogeneity(dat_harm_stg_sugar) #Heterogeneity Q and Q p-val
plt_comb_stg_sugar<-mr_pleiotropy_test(dat_harm_stg_sugar)  #Horizontal pleiotropy
sin_comb_stg_sugar<-mr_singlesnp(dat_harm_stg_sugar, all_method=c("mr_ivw", "mr_egger_regression"))
p_stg_sugar <- mr_forest_plot(sin_comb_stg_sugar)
p_stg_sugar[[1]]
dir_stg_sugar <- directionality_test(dat_harm_stg_sugar)
leave_comb_stg_sugar<- mr_leaveoneout(dat_harm_stg_sugar, parameters = default_parameters())
p_leave_stg_sugar <- mr_leaveoneout_plot(leave_comb_stg_sugar)
p_leave_stg_sugar[[1]]
#
mrpresso_stg_sugar<- run_mr_presso(dat_harm_stg_sugar, NbDistribution = 10000, SignifThreshold = 0.05)
saveRDS(mrpresso_stg_sugar, file = "YOUR_PATHWAY_DIRECTORY/mrpresso_stg_sugar.rds")
###
#save.image(file = "mr_macronutrients.RData")
############################ END ##############################
