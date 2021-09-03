########################################################################
# R code for performing Mendelian randomization 
#######################################################################
# MR of macrnutrient on T2D #
#######################################################################
library(MendelianRandomization)
library(TwoSampleMR)
library(rio)
library(tidyverse)

#Read exposures
#Proteins
prot_fil <- import("Diet_Protein_GWAS_MA_SSGAC_2020_MolPsych.txt")
# format
prot_exp_dat <- format_data(prot_fil_names, type="exposure")
head(prot_exp_dat)
# clump
prot_exp_dat <- clump_data(prot_exp_dat)

#FAT
fat_fil <- import("Diet_Fat_GWAS_MA_SSGAC_2020_MolPsych.txt")
# format
carb_exp_dat <- format_data(carb_fil_names, type="exposure")
head(carb_exp_dat)
# clump
carb_exp_dat <- clump_data(carb_exp_dat)

#SUGAR
sugar_fil <- import("Diet_Sugar_GWAS_MA_SSGAC_2020_MolPsych.txt")
# format
sugar_exp_dat <- format_data(sugar_fil_names, type="exposure")
head(sugar_exp_dat)
# clump
sugar_exp_dat <- clump_data(sugar_exp_dat)
#################################################################################################
## example for protein intake and CHD
#################################################################################################
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

dat_harm_chd_prot<- harmonise_data(
    exposure_dat = prot_exp_dat,
    outcome_dat = chd_out_dat_prot)

export(dat_harm_chd_prot, "Results/dat_harm_chd_prot.xlsx")

# create input
MRInput_chd_prot <- mr_input(bx = dat_harm_chd_prot$beta.exposure,
                             bxse = dat_harm_chd_prot$se.exposure,
                             by = dat_harm_chd_prot$beta.outcome,
                             byse = dat_harm_chd_prot$se.outcome)
#Run MR
MRAll_chd_prot <- mr_allmethods(MRInput_chd_prot, method = "all")

# RAPS
raps_dat_harm_chd_prot <- mr(dat_harm_chd_prot,  method_list = c("mr_raps"))

#Plot
p_chd_prot<- mr_plot(mr_allmethods(MRInput_chd_prot, method = "main"))
p_chd_prot +
    xlab("Protein intake (E%)") +
    ylab("Coronary heart disease (CHD)")+
    ggsave(file = "Results/MRAll_chd_prot_plot.pdf")

# Sensitivity analyses
het_comb_chd_prot<-mr_heterogeneity(dat_harm_chd_prot) #Heterogeneity Q and Q p-val
plt_comb_chd_prot<-mr_pleiotropy_test(dat_harm_chd_prot)  #Horizontal pleiotropy
sin_comb_chd_prot<-mr_singlesnp(dat_harm_chd_prot, all_method=c("mr_ivw", "mr_egger_regression"))
dir_chd_prot <- directionality_test(dat_harm_chd_prot)
leave_comb_chd_prot<- mr_leaveoneout(dat_harm_chd_prot, parameters = default_parameters())

# PRESSO
library(MRPRESSO)
mrpresso_chd_prot<- run_mr_presso(dat_harm_chd_prot, NbDistribution = 10000, SignifThreshold = 0.05)
saveRDS(mrpresso_chd_prot, file = "mr_presso/mrpresso_chd_prot.rds")

# Calculate pseudo R^2 and F-statistic
###########################################################
# calculate 'r.exposure'
lor <- dat_harm_chd_prot$beta.exposure
af <- dat_harm_chd_prot$eaf.exposure
ncase <- dat_harm_chd_prot$ncase.exposure
ncontrol <- dat_harm_chd_prot$ncontrol.exposure
prevalence <- rep(c(0.10), times = nrow(dat_harm_chd_prot)) # prevalence of CHD from literature

r.exposure <- get_r_from_lor(lor, af, ncase, ncontrol, prevalence, model = "logit",
                             correction = FALSE) 

# calculate r^2 and F-statistic
r2_exp <- sum(r.exposure^2)
# f_stat <- (r2_exp*(n-1-k))/((1-r2_exp)*k)  # where n = samplesize.exposure and k = number of IVs
f_stat <- (r2_exp*(12225-1-nrow(dat_harm_chd_prot)))/((1-r2_exp)*nrow(dat_harm_chd_prot))
r2_exp # R^2 for exposure 
f_stat # F-statistic 
###################### END ###########################################
