## Bulk tissue pulse labeling analysis script
source('~/Desktop/Projects/Miceotopes/Bulk/inVivo_func.R')


# Read in searched MS data

#Trypsin_path = "/Users/andrewleduc/Library/CloudStorage/GoogleDrive-research@slavovlab.net/My Drive/MS/Users/aleduc/miceotopes/Bulk/Trypsin/all_search.tsv"
Trypsin_path = 'https://zenodo.org/records/14827610/files/report_trypsin.tsv?download=1'
#LysC_path = "/Users/andrewleduc/Library/CloudStorage/GoogleDrive-research@slavovlab.net/My Drive/MS/Users/aleduc/miceotopes/Bulk/LysC/report.tsv"
LysC_path = 'https://zenodo.org/records/14827610/files/report_lysc.tsv?download=1'

columns_to_read <-c('Genes','Run','Lib.PG.Q.Value','RT','Precursor.Id','Stripped.Sequence','Precursor.Mz',
                    'Precursor.Charge','Precursor.Quantity','Ms1.Area','Protein.Group','Translated.Q.Value','Channel.Q.Value')

Bulk_trypsin <- data.table::fread(Trypsin_path,select = columns_to_read)
Bulk_trypsin <- Bulk_trypsin %>% filter(Lib.PG.Q.Value < .01)
Bulk_trypsin <- Bulk_trypsin %>% filter(Translated.Q.Value < .01)

Bulk_LysC <- data.table::fread(LysC_path,select = columns_to_read)
Bulk_LysC <- Bulk_LysC %>% filter(Lib.PG.Q.Value < .01)
Bulk_LysC <- Bulk_LysC %>% filter(Translated.Q.Value < .01)

# Meta data
#anno = read.csv('/Users/andrewleduc/Desktop/Projects/Miceotopes/Bulk/link.csv')
anno = read.csv('https://zenodo.org/records/14827610/files/Miceotopes_file_annotations.csv?download=1')

####################################################################################
### Make small fasta to reduce search space for recycling (silac variable mod) search
####################################################################################
# library(Biostrings)
# # Read the FASTA file
# fasta_sequences <- readAAStringSet("/Users/andrewleduc/Desktop/Github/QuantQC/inst/extdata/Mouse.fasta")
# 
# # Define the list of protein names to keep
# protein_names <- unique(Bulk_trypsin$Protein.Group)
# 
# # Filter the sequences
# names(fasta_sequences) <- str_extract(names(fasta_sequences), "(?<=\\|)[^|]+(?=\\|)")
# filtered_sequences <- fasta_sequences[names(fasta_sequences) %in% protein_names]
# 
# # Write the filtered sequences to a new FASTA file
# writeXStringSet(filtered_sequences, filepath = "/Users/andrewleduc/Library/CloudStorage/GoogleDrive-research@slavovlab.net/.shortcut-targets-by-id/1uQ4exoKlaZAGnOG1iCJPzYN3ooYYZB7g/MS/Users (1)/aleduc/miceotopes/SC_variable_mod/sc_recycle_version.fasta")


####################################################################################
### Compute AA exchange dynamics from silac variable mod search
####################################################################################

# Data analysis recycling
columns_to_read <-c('Genes','Run','Lib.PG.Q.Value','RT','Precursor.Id','Stripped.Sequence',
                    'Precursor.Charge','Precursor.Quantity','Ms1.Area','Protein.Group','Translated.Q.Value')

# Data from variable mod search
#Recycle <- data.table::fread('/Users/andrewleduc/Library/CloudStorage/GoogleDrive-research@slavovlab.net/My Drive/MS/Users/aleduc/miceotopes/Bulk/Trypsin/searched/report_recycle.tsv')
Recycle <- data.table::fread('https://zenodo.org/records/14827610/files/report_recycle.tsv?download=1')
avg_recycle_estimates <- VarLysinePreproc(Recycle)


avg_recycle_estimates <- avg_recycle_estimates %>% left_join(anno, by = c('Run'))

avg_recycle_estimates$Ratio[avg_recycle_estimates$Tissue == 'BM'] <- avg_recycle_estimates$Ratio[avg_recycle_estimates$Tissue == 'BM'] +.07

# Solve for amino acid exchange from light -> heavy as a func t... H_aa(t)
recycle_df_final <- fit_sum_exp_model(avg_recycle_estimates)

look <- avg_recycle_estimates %>% filter(Time == 5)
ggplot(look,aes(x = Tissue,y=Ratio)) + geom_point(size = 4)



# Plot some of the fits to make sure it looks normal
t <- seq(0,20,by = .1)
n = 6

recycle_df_final$tissue_save[6]
t=c(5,10)
L_aa_fitted <- recycle_df_final$c_save[n] * exp(-recycle_df_final$a_save[n] * t) + recycle_df_final$c_save[n] * exp(-recycle_df_final$b_save[n] * t)
log2(-log(L_aa_fitted))

plot(t,L_aa_fitted,cex.axis = 1.5)

View(recycle_df_final)
####################################################################################
### Pre-processing of Trypsin bulk data
####################################################################################

bad <- c('dAL_LysC_L1','dAL_LysC_L10','dAL_LysC_L5','dAL_AT_Long_M2','dAL_AT_Long_M1')
Bulk_trypsin <- Bulk_trypsin %>% filter(!Run %in% bad)



# Pulse was only with Lysine, so separate R and K peptides
Bulk_trypsin_preproc <- Bulk_preproc(Bulk_trypsin)


Bulk_trypsin_proc <- Bulk_halflife_compute(Bulk_trypsin_preproc,anno,recycle_df_final)


#write.csv(log(2)/Bulk_trypsin_proc$HL,'~/Desktop/Projects/Dilution_paper/Trypsin_Kr.csv')
#write.csv(Bulk_trypsin_proc$prot,'~/Desktop/Projects/Dilution_paper/Trypsin_Abundance.csv')



####################################################################################
### Pre-processing of LysC bulk data
####################################################################################

# Correcting accidental swap when running samples
Bulk_LysC$Run[Bulk_LysC$Run == 'dAL_LysC_M1'] <- 'dAL_LysC_B10_'
Bulk_LysC$Run[Bulk_LysC$Run == 'dAL_LysC_B10'] <- 'dAL_LysC_M1'
Bulk_LysC$Run[Bulk_LysC$Run == 'dAL_LysC_B10_'] <- 'dAL_LysC_B10'

# Removing a couple failed sample prep runs
bad <- c('dAL_LysC_L1','dAL_LysC_L10','dAL_LysC_L5','dAL_LysC_M1','dAL_LysC_M2')
Bulk_LysC <- Bulk_LysC %>% filter(!Run %in% bad)



Bulk_LysC_preproc <- Bulk_preproc(Bulk_LysC)
Bulk_LysC_proc <- Bulk_halflife(Bulk_LysC_preproc,anno,recycle_df_final)



write.csv(log(2)/Bulk_LysC_proc$HL,'~/Desktop/Projects/Dilution_paper/LysC_Kr.csv')
write.csv(Bulk_LysC_proc$prot,'~/Desktop/Projects/Dilution_paper/LysC_Abundance.csv')


####################################################################################
### Supplemental figure 2, plots a protein numbers
####################################################################################

# Plot coverage stats 
count_stats(Bulk_trypsin_proc,bulk_total_prot_trypsin,Bulk_LysC_proc)
plot_precursors(Bulk_LysC,Bulk_trypsin)



####################################################################################
### Supplemental figure 2, plots b (Abundance accuracy lysc trypsin digest comp)
####################################################################################

sect <- intersect(rownames(Bulk_LysC_proc$prot),rownames(Bulk_trypsin_proc$prot))

# adding some updated information to annotations file
anno$dig <- NA
anno$dig[anno$Run %in% colnames(Bulk_LysC_proc$prot)] <- 'LysC'
anno$dig[anno$Run %in% colnames(Bulk_trypsin_proc$prot)] <- 'Trypsin'
anno <- anno %>% filter(Run %in% c(colnames(Bulk_LysC_proc$prot),colnames(Bulk_trypsin_proc$prot)))
anno$Sample <- str_extract(anno$Run, "(?<=_)[^_]+$") 


pf <- rowSums(is.na(Bulk_trypsin_proc$prot[sect,anno$Run[anno$dig == 'Trypsin' & anno$Tissue == 'BM']]))
hist(pf)
pf <- pf[pf ==0]
sect <- intersect(sect,names(pf))

df_abs_brain <- data.frame(Trypsin = log2(rowMeans(Bulk_trypsin_proc$prot[sect,anno$Run[anno$dig == 'Trypsin' & anno$Tissue == 'Brain']],na.rm=T)),
                           LysC = log2(rowMeans(Bulk_LysC_proc$prot[sect,anno$Run[anno$dig == 'LysC' & anno$Tissue == 'Brain']],na.rm=T)))
df_abs_brain$Tissue = 'Brain'
df_abs_lung <- data.frame(Trypsin = log2(rowMeans(Bulk_trypsin_proc$prot[sect,anno$Run[anno$dig == 'Trypsin' & anno$Tissue == 'Lung']],na.rm=T)),
                           LysC = log2(rowMeans(Bulk_LysC_proc$prot[sect,anno$Run[anno$dig == 'LysC' & anno$Tissue == 'Lung']],na.rm=T)))
df_abs_lung$Tissue = 'Lung'
df_abs_BM <- data.frame(Trypsin = log2(rowMeans(Bulk_trypsin_proc$prot[sect,anno$Run[anno$dig == 'Trypsin' & anno$Tissue == 'BM']],na.rm=T)),
                          LysC = log2(rowMeans(Bulk_LysC_proc$prot[sect,anno$Run[anno$dig == 'LysC' & anno$Tissue == 'BM']],na.rm=T)))
df_abs_BM$Tissue = 'BM'

df_abs <- rbind(df_abs_BM,df_abs_lung,df_abs_brain)

ggplot(df_abs,aes(x = 2^LysC,y=2^Trypsin)) + geom_point(alpha = .2)+theme_bw()+
  facet_wrap(~Tissue)+ scale_y_log10()+ scale_x_log10()



####################################################################################
### Supplemental figure 2, plots c (Clearance rate, accuracy lysc trypsin digest comp)
####################################################################################

# Find common space of proteins with deg rate quantified in each
sect <- intersect(rownames(Bulk_LysC_proc$HL),rownames(Bulk_trypsin_proc$HL))

# Grab data out of processed data lists
halflife_lysC_prot <- Bulk_LysC_proc$HL
halflife_df_prot_trypsin <- Bulk_trypsin_proc$HL
halflife_lysC_prot2 <- halflife_lysC_prot
halflife_df_prot_trypsin2 <- halflife_df_prot_trypsin

bulk_total_lysc_prot <- Bulk_LysC_proc$prot
bulk_total_prot_trypsin <- Bulk_trypsin_proc$prot
bulk_total_prot_trypsin2 <- bulk_total_prot_trypsin
bulk_total_lysc_prot2 <- bulk_total_lysc_prot

# Create common sample identifiers
colnames(halflife_lysC_prot) <- str_split(colnames(halflife_lysC_prot), "_", simplify = TRUE)[, 3] 
colnames(halflife_df_prot_trypsin) <- str_split(colnames(halflife_df_prot_trypsin), "_", simplify = TRUE)[, 4]
colnames(bulk_total_prot_trypsin2) <- str_split(colnames(bulk_total_prot_trypsin2), "_", simplify = TRUE)[, 4]
colnames(bulk_total_lysc_prot2) <- str_split(colnames(bulk_total_lysc_prot2), "_", simplify = TRUE)[, 3]


for(i in 1:ncol(halflife_df_prot_trypsin2)){
  halflife_df_prot_trypsin2[,i] <- halflife_df_prot_trypsin2[,i] /median(halflife_df_prot_trypsin2[,i],na.rm = T)
  bulk_total_prot_trypsin2[,i] <- bulk_total_prot_trypsin2[,i] /median(bulk_total_prot_trypsin2[,i],na.rm = T)
  
}


for(i in 1:ncol(halflife_lysC_prot2)){
  halflife_lysC_prot2[,i] <- halflife_lysC_prot2[,i] /median(halflife_lysC_prot2[,i],na.rm = T)
  bulk_total_lysc_prot2[,i] <- bulk_total_lysc_prot2[,i] /median(bulk_total_lysc_prot2[,i],na.rm = T)
}


# Brain

Brain_filt_prot <- Compare_HL_quant(halflife_lysC_prot2,bulk_total_prot_trypsin2, anno,'Brain')
plot(rowMeans(log2(halflife_df_prot_trypsin[sect,Brain_filt_prot$samp]),na.rm = T),rowMeans(log2(halflife_lysC_prot[sect,Brain_filt_prot$samp]),na.rm = T))
plot(rowMeans(log2(halflife_df_prot_trypsin[Brain_filt_prot$prots,Brain_filt_prot$samp]),na.rm = T),rowMeans(log2(halflife_lysC_prot[Brain_filt_prot$prots,Brain_filt_prot$samp]),na.rm = T))
Brain_HL_plot <- data.frame(
  Trypsin = rowMeans(log2(halflife_df_prot_trypsin[Brain_filt_prot$prots, Brain_filt_prot$samp]), na.rm = TRUE),
  LysC = rowMeans(log2(halflife_lysC_prot[Brain_filt_prot$prots, Brain_filt_prot$samp]), na.rm = TRUE)
)
Brain_HL_plot$Tissue <- 'Brain'


# Store SDs
Brain_kr <-sd(rowMeans(log2(halflife_df_prot_trypsin[Brain_filt_prot$prots,Brain_filt_prot$samp]),na.rm = T),na.rm=T)
Brain_kd <- sd(log2(2^rowMeans(log2(halflife_df_prot_trypsin[Brain_filt_prot$prots,Brain_filt_prot$samp]),na.rm = T) - log(2)/median(store$store[store$Tissue == 'Brain'])),na.rm=T)
Brain_kg <- log(2)/median(store$store[store$Tissue == 'Brain'])


# Bone Marrow

BM_filt_prot <- Compare_HL_quant(halflife_lysC_prot2,bulk_total_prot_trypsin2,anno, 'BM')
plot(rowMeans(log2(halflife_df_prot_trypsin[sect,BM_filt_prot$samp]),na.rm = T),rowMeans(log2(halflife_lysC_prot[sect,BM_filt_prot$samp]),na.rm = T))
plot(rowMeans(log2(halflife_df_prot_trypsin[BM_filt_prot$prots,BM_filt_prot$samp]),na.rm = T),rowMeans(log2(halflife_lysC_prot[BM_filt_prot$prots,BM_filt_prot$samp]),na.rm = T))
BM_HL_plot <- data.frame(
  Trypsin = rowMeans(log2(halflife_df_prot_trypsin[BM_filt_prot$prots, BM_filt_prot$samp]), na.rm = TRUE),
  LysC = rowMeans(log2(halflife_lysC_prot[BM_filt_prot$prots, BM_filt_prot$samp]), na.rm = TRUE)
)
BM_HL_plot$Tissue <- 'Bone marrow'

# Store SDs
BM_kr <- sd(rowMeans(log2(halflife_lysC_prot[BM_filt_prot$prots,BM_filt_prot$samp]),na.rm = T),na.rm=T)
BM_kd <- sd(log2(2^rowMeans(log2(halflife_df_prot_trypsin[BM_filt_prot$prots,BM_filt_prot$samp]),na.rm = T) - log(2)/median(store$store[store$Tissue == 'BM'])),na.rm=T)
BM_kg <- log(2)/median(store$store[store$Tissue == 'BM'])


# Lung

Lung_filt_prot <- Compare_HL_quant(halflife_lysC_prot2,bulk_total_prot_trypsin2,anno, 'Lung')
plot(rowMeans(log2(halflife_df_prot_trypsin[sect,Lung_filt_prot$samp]),na.rm = T),rowMeans(log2(halflife_lysC_prot[sect,Lung_filt_prot$samp]),na.rm = T))
plot(rowMeans(log2(halflife_df_prot_trypsin[Lung_filt_prot$prots ,Lung_filt_prot$samp]),na.rm = T),rowMeans(log2(halflife_lysC_prot[Lung_filt_prot$prots ,Lung_filt_prot$samp]),na.rm = T))
Lung_HL_plot <- data.frame(
  Trypsin = rowMeans(log2(halflife_df_prot_trypsin[Lung_filt_prot$prots, Lung_filt_prot$samp]), na.rm = TRUE),
  LysC = rowMeans(log2(halflife_lysC_prot[Lung_filt_prot$prots, Lung_filt_prot$samp]), na.rm = TRUE)
)
Lung_HL_plot$Tissue <- 'Alveolus'

# Store SDs
Lung_kr <- sd(rowMeans(log2(halflife_df_prot_trypsin[Lung_filt_prot$prots ,Lung_filt_prot$samp]),na.rm = T),na.rm=T)
Lung_kd <- sd(log2(2^rowMeans(log2(halflife_df_prot_trypsin[Lung_filt_prot$prots ,Lung_filt_prot$samp]),na.rm = T) - log(2)/median(store$store[store$Tissue == 'Lung'])),na.rm=T)
Lung_kg <- log(2)/median(store$store[store$Tissue == 'Lung'])

# Pancreas... no LysC digest but same filtering

anno_hold <- anno %>% filter(Tissue == 'Pancreas')
halflife_df_prot_trypsin2_hold <- halflife_df_prot_trypsin2[,colnames(halflife_df_prot_trypsin2) %in% anno_hold$Sample]
vect <- rownames(halflife_df_prot_trypsin2) 
vect[rowSums(is.na(halflife_df_prot_trypsin2[,sect_sample])==F) < 4] <- NA
vect[rowSds(halflife_df_prot_trypsin2[,sect_sample],na.rm = T)/rowMeans(halflife_df_prot_trypsin2[,sect_sample],na.rm = T) > .4] <- NA
vect <- vect[is.na(vect) ==F]
Pancreas_filt_prot <- vect

# Store SDs
Panc_kd <- sd(log2(rowMeans(halflife_df_prot_trypsin2[Pancreas_filt_prot,sect_sample],na.rm = T) -  log(2)/median(store$store[store$Tissue == 'Brain'])),na.rm=T)
Panc_kr <- sd(log2(rowMeans(halflife_df_prot_trypsin2[Pancreas_filt_prot,sect_sample],na.rm = T)),na.rm=T)
Panc_kg <- log(2)/median(store$store[store$Tissue == 'Brain'])

plot(log2(rowMeans(halflife_df_prot_trypsin2_hold[Pancreas_filt_prot,c('P1','P2','P3')],na.rm=T)),log2(rowMeans(halflife_df_prot_trypsin2_hold[Pancreas_filt_prot,c('P4','P5','P9')],na.rm=T)))

Panc_HL_plot <- data.frame(
  Trypsin = log2(rowMeans(halflife_df_prot_trypsin2_hold[Pancreas_filt_prot,c('P1','P2','P3')],na.rm=T)),
  LysC = log2(rowMeans(halflife_df_prot_trypsin2_hold[Pancreas_filt_prot,c('P4','P5','P9')],na.rm=T))
)
Panc_HL_plot$Tissue <- 'Pancreas'


##### ##### ##### ##### 
##### This is figure 3C
##### ##### ##### ##### 
df_rates <- data.frame(tissue = c('Brain','Lung','BM','Pancreas'),
                       Kr = c(Brain_kr, Lung_kr,BM_kr,Panc_kr),
                       Kd = c(Brain_kd, Lung_kd,BM_kd,Panc_kd),
                       Kg = c(Brain_kg, Lung_kg,BM_kg,Panc_kg))
df_rates <- reshape2::melt(df_rates, id = c('Kg','tissue'))

ggplot(df_rates, aes(x = Kg,y = value,color = variable)) + geom_point()+
  theme_classic(base_size = 15)
ggplot(df_rates, aes(x = Kd, y = Kr, size = Kg, color = tissue)) + 
  geom_point() +
  theme_classic(base_size = 15) + 
  ylim(c(0.6, 1.1)) + 
  xlim(c(0.75, 1.25)) +
  scale_size_continuous(range = c(4, 15)) +  # Increase point sizes
  theme(
    panel.border = element_rect(color = "black", size = 1.5, fill = NA)  # Add border only around main plot
  )+scale_color_manual(values = c('#BC3EFA','#088109','#FB8A02','#0237FB'))



##### ##### ##### ##### 
##### This is Supplemental fig 2C
##### ##### ##### ##### 
reliabile_deg <- rbind(Lung_HL_plot,BM_HL_plot,Brain_HL_plot,Panc_HL_plot)
ggplot(reliabile_deg,aes(x = log(2)/2^Trypsin,y = log(2)/2^LysC)) + geom_point(alpha = .1) +
  facet_wrap(~Tissue,scales = 'free')+theme_bw() + scale_x_log10()+ scale_y_log10()



# Save all the well quantified clearance rate proteins
Prot_good_List <- list(BM = BM_filt_prot$prots ,Brain = Brain_filt_prot$prots,Lung = Lung_filt_prot$prots,Pancreas =Pancreas_filt_prot )



####################################################################################
### Figure 3, D-F PCA analysis of clearance rates
####################################################################################


#### #### #### #### #### #### #### #### #### #### #### 
#### Bone marrow vs brain comparison
#### #### #### #### #### #### #### #### #### #### #### 

secttt <- intersect(BM_filt_prot$prots, Brain_filt_prot$prots)

# take average Kr from both digests
M1 <- rowMeans(cbind(rowMeans(log2(halflife_df_prot_trypsin[secttt ,BM_filt_prot$samp]),na.rm = T),rowMeans(log2(halflife_lysC_prot[secttt ,BM_filt_prot$samp]),na.rm = T)),na.rm=T)
B1 <- rowMeans(cbind(rowMeans(log2(halflife_df_prot_trypsin[secttt ,Brain_filt_prot$samp]),na.rm = T),rowMeans(log2(halflife_lysC_prot[secttt ,Brain_filt_prot$samp]),na.rm = T)),na.rm=T)
M1 <- log2(log(2)/2^M1)
B1 <- log2(log(2)/2^B1)

# What is the reliability of the Lung - brain Kr fold change
MovB1_kr <- rowMeans(log2(halflife_df_prot_trypsin[secttt ,BM_filt_prot$samp]),na.rm = T) - rowMeans(log2(halflife_df_prot_trypsin[secttt ,Brain_filt_prot$samp]),na.rm = T)
MovB2_kr <- rowMeans(log2(halflife_lysC_prot[secttt ,BM_filt_prot$samp]),na.rm = T) - rowMeans(log2(halflife_lysC_prot[secttt ,Brain_filt_prot$samp]),na.rm = T)
cor(MovB1_kr,MovB2_kr)

# take average concentration from both digests
M1_p <- rowMeans(cbind(rowMeans(log2(bulk_total_prot_trypsin2[secttt ,BM_filt_prot$samp]),na.rm = T),rowMeans(log2(bulk_total_lysc_prot2[secttt ,BM_filt_prot$samp]),na.rm = T)),na.rm=T)
B1_p <- rowMeans(cbind(rowMeans(log2(bulk_total_prot_trypsin2[secttt ,Brain_filt_prot$samp]),na.rm = T),rowMeans(log2(bulk_total_lysc_prot2[secttt ,Brain_filt_prot$samp]),na.rm = T)),na.rm=T)

# What is the reliability of the Lung - brain concentration fold change
MovB1_p <- rowMeans(log2(bulk_total_prot_trypsin2[secttt ,BM_filt_prot$samp]),na.rm = T) - rowMeans(log2(bulk_total_prot_trypsin2[secttt ,Brain_filt_prot$samp]),na.rm = T)
MovB2_p <- rowMeans(log2(bulk_total_lysc_prot2[secttt ,BM_filt_prot$samp]),na.rm = T) - rowMeans(log2(bulk_total_lysc_prot2[secttt ,Brain_filt_prot$samp]),na.rm = T)
cor(MovB1_p,MovB2_p,use = 'pairwise.complete.obs')

reliability_Brain_v_marrow <- cor(MovB1_p,MovB2_p,use = 'pairwise.complete.obs')*cor(MovB1_kr,MovB2_kr)

## Figure 3 D, part1
df.shared <- data.frame(BM = (M1) , Brain =  (B1))
ggplot(df.shared,aes(x = Brain,y = BM)) + geom_point(size = 3,alpha = .2) +
  dot_plot  + geom_abline(slope = 1/TLS(B1,M1)[[1]], intercept = -.5) +
  xlim(c(-6,-1)) + ylim(c(-6,1)) + ylab('Bone marrow')


# PCA analysis
data_matrix <- cbind(B1-mean(B1),M1-mean(M1))
pca_result <- prcomp(data_matrix, center = TRUE, scale. = FALSE)

# Get standard deviations (singular values)
std_devs <- pca_result$sdev

# Compute variance explained
var_explained <- (std_devs^2) / sum(std_devs^2) * 100

## Figure 3 D, part2
pca_result <- as.data.frame(pca_result$x)
ggplot(pca_result,aes(x=PC1,y=PC2))+
  geom_point(size = 3,alpha = .2) +
  dot_plot+ ylim(c(-3,3)) + xlim(c(-3,3))


## Figure 3 F part 1 
df.shared <- data.frame(PC1 = pca_result$PC1,PC2 = pca_result$PC2, Diff =  M1 - B1,abs_dif = M1_p - B1_p)
ggplot(df.shared,aes(x=PC2,y=Diff))+
  geom_point(size = 3,alpha = .2) +
  dot_plot
cor(df.shared$PC2,df.shared$Diff)^2

## Figure 3 F part 2 
ggplot(df.shared,aes(x=PC1,y=abs_dif))+
  geom_point(size = 3,alpha = .2) +
  dot_plot + xlim(c(-2,2))
cor(df.shared$PC1,df.shared$abs_dif)^2/reliability_Brain_v_marrow



#### #### #### #### #### #### #### #### #### #### #### 
#### Lung vs brain comparison
#### #### #### #### #### #### #### #### #### #### #### 
secttt <- intersect(Lung_filt_prot$prots, Brain_filt_prot$prots)

# take average Kr from both digests
L1 <- rowMeans(cbind(rowMeans(log2(halflife_df_prot_trypsin[secttt ,Lung_filt_prot$samp]),na.rm = T),rowMeans(log2(halflife_lysC_prot[secttt ,Lung_filt_prot$samp]),na.rm = T)),na.rm=T)
B1 <- rowMeans(cbind(rowMeans(log2(halflife_df_prot_trypsin[secttt ,Brain_filt_prot$samp]),na.rm = T),rowMeans(log2(halflife_lysC_prot[secttt ,Brain_filt_prot$samp]),na.rm = T)),na.rm=T)
L1 <- log2(log(2)/2^L1)
B1 <- log2(log(2)/2^B1)

# What is the reliability of the Lung - brain Kr fold change
LovB1_kr <- rowMeans(log2(halflife_df_prot_trypsin[secttt ,Lung_filt_prot$samp]),na.rm = T) - rowMeans(log2(halflife_df_prot_trypsin[secttt ,Brain_filt_prot$samp]),na.rm = T)
LovB2_kr <- rowMeans(log2(halflife_lysC_prot[secttt ,Lung_filt_prot$samp]),na.rm = T) - rowMeans(log2(halflife_lysC_prot[secttt ,Brain_filt_prot$samp]),na.rm = T)
cor(LovB1,LovB2)

# take average concentration from both digests
L1_p <- rowMeans(cbind(rowMeans(log2(bulk_total_prot_trypsin2[secttt ,Lung_filt_prot$samp]),na.rm = T),rowMeans(log2(bulk_total_lysc_prot2[secttt ,Lung_filt_prot$samp]),na.rm = T)),na.rm=T)
B1_p <- rowMeans(cbind(rowMeans(log2(bulk_total_prot_trypsin2[secttt ,Brain_filt_prot$samp]),na.rm = T),rowMeans(log2(bulk_total_lysc_prot2[secttt ,Brain_filt_prot$samp]),na.rm = T)),na.rm=T)

# What is the reliability of the Lung - brain concentration fold change
LovB1_p <- rowMeans(log2(bulk_total_prot_trypsin2[secttt ,Lung_filt_prot$samp]),na.rm = T) - rowMeans(log2(bulk_total_prot_trypsin2[secttt ,Brain_filt_prot$samp]),na.rm = T)
LovB2_p <- rowMeans(log2(bulk_total_lysc_prot2[secttt ,Lung_filt_prot$samp]),na.rm = T) - rowMeans(log2(bulk_total_lysc_prot2[secttt ,Brain_filt_prot$samp]),na.rm = T)
cor(LovB1_p,LovB2_p,use = 'pairwise.complete.obs')

reliability_Brain_v_lung <- cor(LovB1_p,LovB2_p,use = 'pairwise.complete.obs')*cor(LovB1,LovB2)


## Figure 3 E, part1
df.shared <- data.frame(BM = (L1) , Brain =  (B1))
ggplot(df.shared,aes(x = Brain,y = BM)) + geom_point(size = 3,alpha = .2) +
  dot_plot  + geom_abline(slope = 1/TLS(B1,L1)[[1]], intercept = 1.1) +
  xlim(c(-6,-1)) + ylim(c(-6,1)) + ylab('Alveolus')


# PCA analysis
data_matrix <- cbind(B1-mean(B1),L1-mean(L1))
pca_result <- prcomp(data_matrix, center = TRUE, scale. = FALSE)

# Get standard deviations (singular values)
std_devs <- pca_result$sdev

# Compute variance explained
var_explained <- (std_devs^2) / sum(std_devs^2) * 100

## Figure 3 E, part2
pca_result <- as.data.frame(pca_result$x)
ggplot(pca_result,aes(x=PC1,y=PC2))+
  geom_point(size = 3,alpha = .2) +
  dot_plot+ ylim(c(-3,3)) + xlim(c(-3,3))


## Figure 3 G part 1 
df.shared <- data.frame(PC1 = pca_result$PC1,PC2 = pca_result$PC2, Diff =  L1 - B1,abs_dif = L1_p - B1_p)
ggplot(df.shared,aes(x=PC2,y=Diff))+
  geom_point(size = 3,alpha = .2) +
  dot_plot
cor(df.shared$PC2,df.shared$Diff)^2

## Figure 3 G part 2 
ggplot(df.shared,aes(x=PC2,y=abs_dif))+
  geom_point(size = 3,alpha = .2) +
  dot_plot + xlim(c(-2,2))
cor(df.shared$PC2,df.shared$abs_dif)^2/(0.7113127*0.8302321)






####################################################################################
### Supplemental Figure 2, d Getting cell doubling times (division rates)
####################################################################################

# Get division rate through histone turnover
store <- colMedians(Bulk_trypsin_proc$HL[c('P02301;P84244','P62806','Q64478;Q6ZWY9'),],na.rm=T)
store <- as.data.frame(store)
store$Run <- rownames(store)
store <- store %>% left_join(anno, by = c('Run'))
rownames(store) <- store$Sample

store2 <- colMedians(Bulk_LysC_proc$HL[c('P02301;P84244','P62806','Q64478;Q6ZWY9'),],na.rm=T)
store2 <- as.data.frame(store2)
store2$Run <- rownames(store2)
store2 <- store2 %>% left_join(anno, by = c('Run'))
rownames(store2) <- store2$Sample

store_sect <- intersect(rownames(store),rownames(store2))
plot(store2$store2[store2$Sample %in% store_sect],store$store[store$Sample %in% store_sect])
abline(a=0,b=1)


store <- store %>% left_join(store2, by=c('Sample','Tissue','Age','Gender','Time'))
store$Average <-  rowSums(cbind(store$store, store$store2), na.rm = TRUE)

ggplot(store, aes(x = store,y = store2,color = Tissue))+geom_point(size = 5,alpha = .75) +
  theme_classic(base_size = 15) + scale_x_log10()+scale_y_log10()+
  theme(
    axis.text = element_text(size = 14, color = "black"),  # Larger axis text
    axis.title = element_text(size = 16,),  # Bold axis titles
    legend.title = element_text(size = 14),  # Larger legend title
    legend.text = element_text(size = 12),  # Larger legend text
    panel.border = element_rect(color = "black", fill = NA, size = 1),  # Add top and left border lines
    axis.line.x.top = element_line(color = "black"),  # Top x-axis line
    axis.line.y.left = element_line(color = "black")  # Left y-axis lin
  ) + scale_color_manual(values = c('#BC3EFA','#088109','#FB8A02','#0237FB'))+
  geom_abline(slope = 1) + xlab('Trypsin') + ylab('LysC') +
  ggtitle('Population average doubling time, days')


####################################################################################
### Figure 3, b (division rate by R^2)
####################################################################################


Rsq_store <- c()
Predicted_Rsq <- c()
div_Rate <- c()
corrected <- c()
for(i in 1:nrow(store)){
  div_Rate <- c(div_Rate,store$Average[store$Sample == anno$Sample[i]])
  
  tissue_hold <- store$Tissue[i]
  prot_list <- Prot_good_List[[tissue_hold]]
  
  
  sample_hold <- store$Sample[i]
 
  if(is.na(store$store2[i] + store$store[i])==T){
    corrected <- c(corrected,F)
    
    
    if(is.na(store$Run.x[i])==F){
      Rsq_store <- c(Rsq_store,cor(log2(halflife_df_prot_trypsin2[prot_list,sample_hold]),
                                   log2(bulk_total_prot_trypsin2[prot_list,sample_hold]),use = 'pairwise.complete.obs')^2 / .6)
      
      Predicted_Rsq <- c(Predicted_Rsq,estimate_Rsq(log(2)/halflife_df_prot_trypsin2[prot_list,sample_hold],bulk_total_prot_trypsin2[prot_list,sample_hold]))
    }
    if(is.na(store$Run.y[i])==F){
      
      Rsq_store <- c(Rsq_store,cor(log2(halflife_lysC_prot2[prot_list,sample_hold]),
                                   log2(bulk_total_lysc_prot2[prot_list,sample_hold]),use = 'pairwise.complete.obs')^2)
      Predicted_Rsq <- c(Predicted_Rsq,estimate_Rsq(log(2)/halflife_lysC_prot2[prot_list,sample_hold],bulk_total_lysc_prot2[prot_list,sample_hold]))
    }
    
    
  }else{
    corrected <- c(corrected,T)
    
    HL <- rowMeans(cbind(halflife_lysC_prot2[prot_list,sample_hold],halflife_df_prot_trypsin2[prot_list,sample_hold]),na.rm = T)
    Abs <- rowMeans(cbind(bulk_total_lysc_prot2[prot_list,sample_hold],bulk_total_prot_trypsin2[prot_list,sample_hold]),na.rm = T)
    
    
    
    Rsq <- cor(log2(HL),log2(Abs),use = 'pairwise.complete.obs')
    if(Rsq < .0){
      Rsq <- .01
    }else{
      Rsq <- Rsq^2
      
    }
    
    correction1 <- cor(log2(halflife_lysC_prot2[prot_list,sample_hold]),log2(halflife_df_prot_trypsin2[prot_list,sample_hold]),use = 'pairwise.complete.obs')
    correction2 <- cor(log2(bulk_total_lysc_prot2[prot_list,sample_hold]),log2(bulk_total_prot_trypsin2[prot_list,sample_hold]),use = 'pairwise.complete.obs')
    
    Rsq_store <- c(Rsq_store,Rsq/(correction1*correction2))
    
    Predicted_Rsq <- c(Predicted_Rsq,estimate_Rsq(log(2)/HL,Abs,tissue_hold))
    
    
  }
}



df_final <- data.frame(Rsq_empirical = Rsq_store,Rsq_predicted = Predicted_Rsq,
                       spearman = corrected,div_rate = div_Rate,Tissue = store$Tissue)


df_final <- df_final %>% filter(Rsq_predicted < .3)

# Bone marrow rates were slightly correlated in wrong direction which would anti-contribute 
# so I just set Rsq to a minimum value
df_final$Rsq_empirical[df_final$Tissue == 'BM'] <- .01


# 
ggplot(df_final,aes(x = div_rate,y = Rsq_empirical,color = Tissue))+geom_point(size = 5,alpha = .75) +
  ylim(c(0,.5)) + theme_classic(base_size = 15) +
  theme(
    axis.text = element_text(size = 14, color = "black"),  # Larger axis text
    axis.title = element_text(size = 16,),  # Bold axis titles
    legend.title = element_text(size = 14),  # Larger legend title
    legend.text = element_text(size = 12),  # Larger legend text
    panel.border = element_rect(color = "black", fill = NA, size = 1),  # Add top and left border lines
    axis.line.x.top = element_line(color = "black"),  # Top x-axis line
    axis.line.y.left = element_line(color = "black")  # Left y-axis lin
  ) + scale_color_manual(values = c('#BC3EFA','#088109','#FB8A02','#0237FB')) + 
  scale_x_log10(limits = c(2, 200))



# ggplot(df_final,aes(x = Rsq_empirical,y = Rsq_predicted))+
#   geom_point(size = 5) +ylim(c(0,.3))+theme_classic(base_size = 15) +
#   theme(
#     axis.text = element_text(size = 14, color = "black"),  # Larger axis text
#     axis.title = element_text(size = 16,),  # Bold axis titles
#     legend.title = element_text(size = 14),  # Larger legend title
#     legend.text = element_text(size = 12),  # Larger legend text
#     panel.border = element_rect(color = "black", fill = NA, size = 1),  # Add top and left border lines
#     axis.line.x.top = element_line(color = "black"),  # Top x-axis line
#     axis.line.y.left = element_line(color = "black")  # Left y-axis lin
#   )



####################################################################################
### Figure 4, a half life diff vs prot dif
####################################################################################

# Just plotting the Kr rates without any transformations against each other

anno$dig <- NA
anno$dig[grepl('Lys',anno$Run) == T] <- 'LysC'
anno$dig[grepl('Lys',anno$Run) == F] <- 'Trypsin'
anno_brain <- anno %>% filter(Tissue == 'Brain' & dig == 'Trypsin')
anno_panc <- anno %>% filter(Tissue == 'Pancreas')

store_brain <- store %>% filter(Tissue == 'Brain')
store_panc <- store %>% filter(Tissue == 'Pancreas')

plot_sect <- intersect(Prot_good_List$Brain,Prot_good_List$Pancreas)

brain_kr <- log(2)/rowMeans(halflife_df_prot_trypsin[plot_sect,anno_brain$Sample],na.rm=T)
brain_kd <- brain_kr - log(2)/median(store_brain$store)

panc_kr <- log(2)/rowMeans(halflife_df_prot_trypsin[plot_sect,anno_panc$Sample],na.rm=T)
panc_kd <- panc_kr - log(2)/median(store_panc$store)

kr_dif <- log2(brain_kr/panc_kr) - median(log2(brain_kr/panc_kr),na.rm=T)
kd_dif <- log2(brain_kd/panc_kd) - median(log2(brain_kd/panc_kd),na.rm=T)

plot_kr <- data.frame(brain = log2(brain_kr) - median(log2(brain_kr),na.rm=T),panc = log2(panc_kr) - median(log2(panc_kr),na.rm=T))
ggplot(plot_kr, aes(x = brain , y = panc)) + geom_point() +theme_classic(base_size = 15)+
  geom_abline(intercept = 0,slope = 1)



####################################################################################
### Figure 4, b,c (Relative difference explained by deg R^2)
####################################################################################

# Getting proteins with well quantified half life in each tissue for comparison
sect_good <- intersect(intersect(intersect(Prot_good_List$BM,Prot_good_List$Brain),
                       Prot_good_List$Lung),Prot_good_List$Pancreas)


# For storing half life and protein abundance fold changes
cor_matrix <- matrix(NA, nrow = ncol(halflife_lysC_prot2), ncol = ncol(bulk_total_lysc_prot2))
cor_matrix <- matrix(NA, nrow = ncol(halflife_df_prot_trypsin2), ncol = ncol(halflife_df_prot_trypsin2))
# Populate the matrix
for (i in 1:ncol(halflife_df_prot_trypsin2)) {
  for (j in 1:ncol(halflife_df_prot_trypsin2)) {
    P_dif <- log2(bulk_total_prot_trypsin2[sect_good, i] / bulk_total_prot_trypsin2[sect_good, j])
    HL_dif <- log2(halflife_df_prot_trypsin2[sect_good, i] / halflife_df_prot_trypsin2[sect_good, j])
    
    cor_matrix[i, j] <- cor(HL_dif, P_dif, use = 'pairwise.complete.obs')
  }
}


# making sure sample names are on col/rows
colnames(cor_matrix) <- colnames(halflife_df_prot_trypsin2)
rownames(cor_matrix) <- colnames(halflife_df_prot_trypsin2)

# NA repeated upper tri
cor_matrix[upper.tri(cor_matrix)] <- NA

# Melt it for plotting and make upper try and diagonal go away so only have each
# pairwise comparison once
cor_matrix.m <- melt(cor_matrix)
cor_matrix.m$Var1 <- stringr::str_sub(cor_matrix.m$Var1 ,0,1)
cor_matrix.m$Var2 <- stringr::str_sub(cor_matrix.m$Var2 ,0,1)
cor_matrix.m$Var3 <- paste0(cor_matrix.m$Var1,cor_matrix.m$Var2)
cor_matrix.m <- cor_matrix.m %>% filter(is.na(value)==F)


cor_matrix.m <- cor_matrix.m %>% group_by(Var3) %>% summarise(value = median(value))

cor_matrix.m$Var2 <- str_sub(cor_matrix.m$Var3 ,0,1)
cor_matrix.m$Var3 <- str_sub(cor_matrix.m$Var3 ,2,2)

# This is Figure 4C
ggplot(cor_matrix.m,aes(y = Var2,x = Var3,fill = value))+ geom_point(shape = 21, stroke = 2,size = 15)+
  theme_classic(base_size = 15)+
  theme(
    axis.text = element_text(size = 14, color = "black"),  # Larger axis text
    axis.title = element_text(size = 16,),  # Bold axis titles
    legend.title = element_text(size = 14),  # Larger legend title
    legend.text = element_text(size = 12),  # Larger legend text
    panel.border = element_rect(color = "black", fill = NA, size = 2),  # Add top and left border lines
    axis.line.x.top = element_line(color = "black"),  # Top x-axis line
    axis.line.y.left = element_line(color = "black")  # Left y-axis lin
  ) + scale_fill_gradient2(midpoint = .2,high = 'red',low = 'white')



# This is just the brain pancreas comp for Figure 4B

# Well quantified protein Kr between brain and panc 
ssss <- intersect(Prot_good_List$Pancreas,Prot_good_List$Brain)
P_dif <- log2(rowMeans(bulk_total_prot_trypsin2[ssss, 24:33],na.rm = T) / rowMeans(bulk_total_prot_trypsin2[ssss, 1:10],na.rm=T))
HL_dif <- log2(rowMeans(halflife_df_prot_trypsin2[ssss, 24:33],na.rm=T) / rowMeans(halflife_df_prot_trypsin2[ssss, 1:10],na.rm = T))

df_PHL <- as.data.frame(cbind(P_dif,HL_dif))
ggplot(df_PHL,aes(y = P_dif,x = log2(log(2)/2^HL_dif))) + geom_point(size = 2)+theme_classic(base_size = 15)+
  theme(
    axis.text = element_text(size = 14, color = "black"),  # Larger axis text
    axis.title = element_text(size = 16,),  # Bold axis titles
    legend.title = element_text(size = 14),  # Larger legend title
    legend.text = element_text(size = 12),  # Larger legend text # Add top and left border lines
    axis.line.x.top = element_line(color = "black"),  # Top x-axis line
    axis.line.y.left = element_line(color = "black")  # Left y-axis lin
  )+xlim(c(-7,6)) + ylab('Abundance') + xlab('k_r') + ggtitle('log2(Brain/Pancreas)')



####################################################################################
### Figure 4, d (Individual proteins R^2)
####################################################################################


# Normalizing all the data to relative log2 fold changes for the plotting
halflife_df_prot_trypsin2_norm <- QuantQC::Normalize_reference_vector(halflife_df_prot_trypsin2,log=T)
bulk_total_prot_trypsin2_norm <- QuantQC::Normalize_reference_vector(bulk_total_prot_trypsin2,log=T)
halflife_lysC_prot2[,i] <- halflife_lysC_prot2[,i] /median(halflife_lysC_prot2[,i],na.rm = T)

halflife_lysC_prot2_norm <-  QuantQC::Normalize_reference_vector(halflife_df_prot_trypsin2,log=T)
bulk_total_lysc_prot2_norm <- QuantQC::Normalize_reference_vector(bulk_total_lysc_prot2,log=T)
sect <- intersect(rownames(bulk_total_lysc_prot2_norm),rownames(halflife_lysC_prot2_norm))
halflife_lysC_prot2_norm <- halflife_lysC_prot2_norm[sect,]
bulk_total_lysc_prot2_norm <- bulk_total_lysc_prot2_norm[sect,]
halflife_lysC_prot2_norm <- halflife_lysC_prot2_norm[,intersect(colnames(bulk_total_lysc_prot2_norm),colnames(halflife_lysC_prot2_norm))]
bulk_total_lysc_prot2_norm <- bulk_total_lysc_prot2_norm[,intersect(colnames(bulk_total_lysc_prot2_norm),colnames(halflife_lysC_prot2_norm))]


# To run on only proteins quantified in both lysC and trypsin
sect <- intersect(rownames(halflife_df_prot_trypsin2_norm),rownames(halflife_lysC_prot2_norm))
#sect2 <- intersect(colnames(halflife_df_prot_trypsin2_norm),colnames(halflife_lysC_prot2_norm))

cor_store_trypsin <- c()
cor_store_LysC <- c()
scramble <- c()
prot_save <- c()
for(i in sect){
  
  if(pairwiseCount(halflife_df_prot_trypsin2_norm[i,sect2],bulk_total_prot_trypsin2_norm[i,sect2]) >15){
    
    cor_store_trypsin <- c(cor_store_trypsin,cor(halflife_df_prot_trypsin2_norm[i,],bulk_total_prot_trypsin2_norm[i,],use = 'pairwise.complete.obs'))
    #cor_store_LysC <- c(cor_store_LysC,cor(bulk_total_lysc_prot2_norm[i,sect2],halflife_lysC_prot2_norm[i,sect2],use = 'pairwise.complete.obs'))
    scramble <- c(scramble,coral(halflife_df_prot_trypsin2_norm[i,],bulk_total_prot_trypsin2_norm[i,]))
  }
  
}



df_prot_cor <- data.frame(prot = sect,cor = cor_store_trypsin, null_dist = scramble)

ggplot(df_prot_cor,aes(x = cor)) +geom_histogram()+theme_classic(base_size = 15) +
  geom_vline(xintercept = 0.35) + xlab('Cor(K_R,Abs)')+ylab('# proteins')


# Select a few example correlations to plot
pplot <- c('P58064','Q60932','O55142','Q9CR68')


           
plot_prot_cor(pplot,anno,bulk_total_prot_trypsin2_norm,halflife_df_prot_trypsin2_norm,bulk_total_lysc_prot2_norm,halflife_lysC_prot2_norm)


####################################################################################
### Figure 4, e and f (Go Terms from across tissue analysis)
####################################################################################

# Mouse go terms
#Go <- read.gmt('/Users/andrewleduc/Desktop/Projects/Miceotopes/Bulk/Gene_sets/m2.all.v2023.2.Mm.symbols.gmt')
Go <- read.gmt('https://zenodo.org/records/14827610/files/GO_Mouse.gmt?download=1')

# Get gene names for mice data
Mouse <- Proc_fasta('/Users/andrewleduc/Desktop/Github/QuantQC/inst/extdata/Mouse.fasta')
df_prot_cor <- df_prot_cor %>% left_join(Mouse, by = c('prot' = 'split_prot'))



# Gene set enrichment analysis on correlation distribution
p_vals <- c()
terms <- c()
med_cor <- c()
med_abs <- c()
for(i in unique(Go$term)){
  Go_hold <- Go %>% filter(term == i)
  
  if(length(intersect(Go_hold$gene,df_prot_cor$split_gene)) > 3){
    df_prot_cor_hold <- df_prot_cor %>% filter(split_gene %in% Go_hold$gene)
    
    p_vals <- c(p_vals,t.test(df_prot_cor_hold$cor,df_prot_cor$null_dist)$p.value)
    terms <- c(terms,i)
    med_cor <- c(med_cor,median(df_prot_cor_hold$cor,na.rm=T))
    
    
    abs_vect <- colMeans(bulk_total_prot_trypsin2_norm[df_prot_cor_hold$prot,],na.rm=T)
    abs_vect <- as.data.frame(abs_vect)
    abs_vect$samp <- rownames(abs_vect)
    
    abs_vect$samp <- stringr::str_sub(abs_vect$samp ,0,1)
    
    abs_vect <- abs_vect %>% group_by(samp) %>% 
      dplyr::summarise(abs_vect = median(abs_vect,na.rm=T))
    
    
    med_abs <- c(med_abs,abs_vect$samp[abs_vect$abs_vect == max(abs_vect$abs_vect)])
    
    
    
    
  }
  
}

# Store results , filter for confidence  
df_go <- data.frame(pval =p_vals,GO = terms,cor = med_cor,highest = med_abs )
df_go$Qval <- p.adjust(p_vals,method = 'BH')
df_go <- df_go %>% filter(Qval < .05)
df_go$Rsq <- df_go$cor^2
#df_go <- df_go %>% filter(Rsq > .2)

# Save stuff for supp tables
write.csv(df_go,'~/Desktop/Projects/Dilution_paper/GO_table3.csv')
write.csv(df_prot_cor,'~/Desktop/Projects/Dilution_paper/Correlations.csv')



#### Saving go terms to be plotted

# Pancreas
big_go <- c('BIOCARTA_EIF2_PATHWAY',
'MORI_PLASMA_CELL_UP',
'REACTOME_INTEGRIN_CELL_SURFACE_INTERACTIONS',
'REACTOME_TRANSCRIPTIONAL_ACTIVATION_OF_MITOCHONDRIAL_BIOGENESIS',
'REACTOME_TRAFFICKING_AND_PROCESSING_OF_ENDOSOMAL_TLR',
'REACTOME_EUKARYOTIC_TRANSLATION_ELONGATION',

# BM 
'REACTOME_MITOTIC_TELOPHASE_CYTOKINESIS',
'REACTOME_DAP12_SIGNALING',
'JACKSON_DNMT1_TARGETS_DN',
'LEE_TARGETS_OF_PTCH1_AND_SUFU_UP',
'REACTOME_PCNA_DEPENDENT_LONG_PATCH_BASE_EXCISION_REPAIR',

# Lung all

'BOYLAN_MULTIPLE_MYELOMA_D_CLUSTER_UP',
'PLASARI_TGFB1_SIGNALING_VIA_NFIC_10HR_DN',
'REACTOME_REGULATION_OF_CYTOSKELETAL_REMODELING_AND_CELL_SPREADING_BY_IPP_COMPLEX_COMPONENTS',
'REACTOME_DEGRADATION_OF_THE_EXTRACELLULAR_MATRIX',

# Brain 

'REACTOME_WNT5A_DEPENDENT_INTERNALIZATION_OF_FZD4',
 'REACTOME_GLYOXYLATE_METABOLISM_AND_GLYCINE_DEGRADATION',
 'WP_SEROTONIN_AND_ANXIETY',
 'REACTOME_COMPLEX_I_BIOGENESIS',
 'BIOCARTA_GSK3_PATHWAY',
 'REACTOME_INSULIN_RECEPTOR_RECYCLING',
 'BIOCARTA_AKAP95_PATHWAY')


## FIGURE 4 E
Go_Terms_boxplot_cor(c('REACTOME_RETROGRADE_NEUROTROPHIN_SIGNALLING','REACTOME_CITRIC_ACID_CYCLE_TCA_CYCLE','BIOCARTA_CK1_PATHWAY','REACTOME_WNT5A_DEPENDENT_INTERNALIZATION_OF_FZD4'),df_prot_cor)


# Reordering for the plotting by Rsq
df_go <- df_go[order(-df_go$Rsq),]
df_go <- df_go %>% filter(GO != i)

# Getting only the go terms we wanted to plot
df_go_plt <- df_go %>% filter(GO %in% big_go)

# Final figure
ggplot(df_go_plt, aes(x = Rsq, y = reorder(GO, Rsq),fill = highest)) +  # Reorder GO by Rsq
  geom_bar(stat = 'identity') +
  ggtitle(i) +
  xlim(c(0, 1)) +
  theme_classic()+
  scale_fill_manual(values = c('#088109','#FB8A02','#BC3EFA','#0237FB'))




