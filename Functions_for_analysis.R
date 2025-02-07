# Metabolic labeling functions
library(stats)
library(matrixStats)
library(stringr)
library(diann)
library(QuantQC)
library(MASS)
library(ggplot2)
library(dplyr)

### Small data analysis functions

TLS <- function(vect1,vect2){

  int_x <- mean(vect1,na.rm = T)
  int_y <- mean(vect2,na.rm = T)
  
  vect1 <- vect1-int_x
  vect2 <- vect2-int_y
  
  mat <- as.matrix(cbind(vect1,vect2))
  
  mat <- mat[complete.cases(mat),]
  
  TLS_mat <- svd(mat)$v
  
  slope <- TLS_mat[1,1]/TLS_mat[2,1]
  
  int <- c(int_x,int_y)
  
  return(list(slope,int))
  
} # Computes total least squares regression slope and intercept

coral <-function(vect1,vect2){
  cor(vect1,vect2,use = 'pairwise.complete.obs',method = 'spearman')
} # Just a spearman correlation that ignors NAs

compute_rsq <- function(P,deg,trans){
  
  
  Y <- P - mean(P,na.rm = T)      #  look$RowMeans1
  X1 <- deg - mean(deg,na.rm = T) 
  
  # Compute predicted values with Beta fixed at 1
  Y_hat <- -X1
  
  # Total Sum of Squares (TSS)
  TSS <- sum((Y - mean(Y,na.rm = T))^2,na.rm = T)
  
  # Residual Sum of Squares (RSS)
  RSS <- sum((Y - Y_hat)^2,na.rm = T)
  
  #RSS1 <- sum((Y - Y_hat))
  
  # R-squared
  rsq1 <- (1 - (RSS / TSS))
  
  
  
  
  
  
  X1 <- trans - mean(trans,na.rm=T) 
  
  # Compute predicted values with Beta fixed at 1
  Y_hat <- X1
  
  # Total Sum of Squares (TSS)
  TSS <- sum((Y - mean(Y,na.rm = T))^2,na.rm = T)
  
  # Residual Sum of Squares (RSS)
  RSS <- sum((Y - Y_hat)^2,na.rm = T)
  
  #RSS2 <- sum((Y - Y_hat))
  
  # R-squared
  rsq2 <- 1 - (RSS / TSS)
  
  rsq1/(rsq2+rsq1)
  
  return(rsq1/(rsq2+rsq1))
  
  #return(RSS1/(RSS2+RSS1))
  
} # This is the way we compute fract variance explained by Kr

count_stats <- function(Bulk_trypsin_proc,bulk_total_prot_trypsin,Bulk_LysC_proc){
  
  
  df_pnumb <- as.data.frame(colSums(is.na(bulk_total_prot_trypsin)==F))
  df_pnumb$Run <- rownames(df_pnumb)
  df_pnumb <- df_pnumb %>% left_join(anno, by = 'Run')
  colnames(df_pnumb)[1] <- 'Protein_numb'
  
  # For counting protein numbers only
  df_pnumb_lysc <- as.data.frame(colSums(is.na(Bulk_LysC_proc$prot)==F))
  df_pnumb_lysc$Run <- rownames(df_pnumb_lysc)
  df_pnumb_lysc <- df_pnumb_lysc %>% left_join(anno, by = 'Run')
  colnames(df_pnumb_lysc)[1] <- 'Protein_numb'
  
  df_pnumb_lysc$dig <- 'LysC'
  df_pnumb$dig <- 'Trypsin'
  df_pnumb <- rbind(df_pnumb,df_pnumb_lysc)
  
  
  # For counting protein numbers with Half life measurments only
  df_tnumb <- as.data.frame(colSums(is.na(Bulk_trypsin_proc$HL)==F))
  df_tnumb$Run <- rownames(df_tnumb)
  df_tnumb <- df_tnumb %>% left_join(anno, by = 'Run')
  colnames(df_tnumb)[1] <- 'Protein_numb'
  df_tnumb$dig <- 'Trypsin'
  
  
  df_tnumb_lysc <- as.data.frame(colSums(is.na(Bulk_LysC_proc$HL)==F))
  df_tnumb_lysc$Run <- rownames(df_tnumb_lysc)
  df_tnumb_lysc <- df_tnumb_lysc %>% left_join(anno, by = 'Run')
  colnames(df_tnumb_lysc)[1] <- 'Protein_numb'
  
  df_tnumb_lysc$dig <- 'LysC'
  df_tnumb <- rbind(df_tnumb,df_tnumb_lysc)
  df_tnumb$Type<- 'Degradation rate'
  df_pnumb$Type<- 'Protein level'
  
  
  
  df_stats <- rbind(df_pnumb,df_tnumb)
  
  
  ggplot(df_stats, aes(x = Tissue, y = Protein_numb, color = dig)) + 
    geom_boxplot(position = position_dodge(width = 0.8), outlier.shape = NA) +  # Dodge boxplots
    geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8), 
                size = 2.5, alpha = 0.7) +  # Dodge and jitter points
    ylim(c(0, 7000)) +
    facet_wrap(~Type, scales = "free_x") +  # Free x-scales if needed
    theme_classic(base_size = 16) +  # Set larger base font size
    theme(
      plot.title = element_text(hjust = 0.5, size = 18),  # Centered and bold title
      axis.text.x = element_text(angle = 45, hjust = 1, size = 14, color = "black"),  # Rotated x-axis labels
      axis.text.y = element_text(size = 14, color = "black"),  # Larger y-axis text
      axis.title.y = element_text(size = 16),  # Larger and bold y-axis title
      legend.title = element_text(size = 14, face = "bold"),  # Larger and bold legend title
      legend.text = element_text(size = 12),  # Larger legend text
      strip.text = element_text(size = 15, face = "bold")  # Larger facet labels
    ) +  # Use a colorblind-friendly palette
    labs(
      x = "",  # Remove x-axis title
      y = "# Proteins quantified",  # Improved y-axis label
    )
}

truncated_sd <- function(x, trim = 0.1, na.rm = TRUE) {
  if (na.rm) x <- x[!is.na(x)]
  trimmed_x <- x[x >= quantile(x, trim) & x <= quantile(x, 1 - trim)]
  return(sd(trimmed_x))
}

robust_sd <- function(x, na.rm = TRUE) {
  if (na.rm) x <- x[!is.na(x)]
  mad <- median(abs(x - median(x)))
  return(mad * 1.4826)  # Scale factor to approximate SD
}

estimate_Rsq <- function(K_r,P,reliability = 'Pancreas'){
  
  if(reliability %in% c('Lung','BM')){
    sigma_logAlpha <- truncated_sd(log2(K_r),trim = .2, na.rm=T)
    #sd(log2(K_r),na.rm=T) 
  }else{
    #sd(log2(K_r),na.rm=T)
    sigma_logAlpha <- sd(log2(K_r), na.rm=T) 
    
  }
  mu_logP <- mean(log2(P*K_r),na.rm=T)
  mu_logAlpha <- mean(log2(K_r),na.rm=T)
  sigma_logP <- sd(log2(P*K_r),na.rm=T)
  
  correlation <- cor(log2(P*K_r),log2(K_r),use = 'pairwise.complete.obs')  # Given correlation
  
  # Covariance matrix
  cov_matrix <- matrix(c(sigma_logP^2, correlation * sigma_logP * sigma_logAlpha, 
                         correlation * sigma_logP * sigma_logAlpha, sigma_logAlpha^2), nrow = 2)
  
  
  simulated_data <- mvrnorm(n = 20000, mu = c(mu_logP, mu_logAlpha), Sigma = cov_matrix)
  
  
  deg <- (simulated_data[, 2])
  trans <- (simulated_data[, 1])
  
  P = trans-deg
  
  return(compute_rsq(P,deg,trans))
  #return(cor(P,deg)^2)
}

Compare_HL_quant <- function(halflife_lysC_prot2,bulk_total_prot_trypsin2, anno, tissue){
  
  i = tissue
  anno_hold <- anno %>% filter(Tissue == i)
  halflife_lysC_prot2_hold <- halflife_lysC_prot2[,colnames(halflife_lysC_prot2) %in% anno_hold$Sample]
  halflife_df_prot_trypsin2_hold <- halflife_df_prot_trypsin2[,colnames(halflife_df_prot_trypsin2) %in% anno_hold$Sample]
  
  
  sect_sample <- intersect(colnames(halflife_lysC_prot2_hold),colnames(halflife_df_prot_trypsin2_hold))
  sect_sample
  
  vect <- rownames(halflife_df_prot_trypsin2) 
  vect[rowSums(is.na(halflife_df_prot_trypsin2[,sect_sample])==F) < 2] <- NA
  vect[rowSds(halflife_df_prot_trypsin2[,sect_sample],na.rm = T)/rowMeans(halflife_df_prot_trypsin2[,sect_sample],na.rm = T) > .5] <- NA
  hist(rowSds(halflife_df_prot_trypsin2[,sect_sample],na.rm = T)/rowMeans(halflife_df_prot_trypsin2[,sect_sample],na.rm = T) )
  vect <- vect[is.na(vect) ==F]
  vect <- intersect(vect,sect)
  length(vect) 
  
  vect2 <- rownames(halflife_lysC_prot2) 
  vect2[rowSums(is.na(halflife_lysC_prot2[,sect_sample])==F) < 2] <- NA
  vect2[rowSds(halflife_lysC_prot2[,sect_sample],na.rm = T)/rowMeans(halflife_lysC_prot2[,sect_sample],na.rm = T) > .5] <- NA
  hist(rowSds(halflife_lysC_prot2[,sect_sample],na.rm = T)/rowMeans(halflife_lysC_prot2[,sect_sample],na.rm = T) )
  vect2 <- vect2[is.na(vect2) ==F]
  length(vect2)  
  vect2 <- intersect(vect,vect2)
  length(vect2)  
  
  #plot(rowMeans(log2(halflife_df_prot_trypsin[vect2,sect_sample]),na.rm = T),rowMeans(log2(halflife_lysC_prot[vect2,sect_sample]),na.rm = T))
  #abline(a = 0,b=1)
  
  return(list(prots = vect2,samp = sect_sample))
}



### Functions for getting Kr by accounting for recycling 

# Preproc of the variable mod search to find peptides with one 
# heavy and one light lysine missed cleaved peptides and taking their ratio to
# regress against time to find the recycling rate function \gamma (t), how avalible light lysine
# changes with time
VarLysinePreproc <- function(Recycle_searched){
  
  # Get double lysine peptides
  Recycle1 <- Recycle_searched %>% filter(str_count(Stripped.Sequence, "K") == 2)
  Recycle1$numb <- stringr::str_count(Recycle1$Precursor.Id, "SILAC")
  Recycle1$seqcharge <- paste0(Recycle1$Stripped.Sequence,Recycle1$Precursor.Charge)
  Recycle1 <- Recycle1 %>% filter(Ms1.Area != 0)
  Recycle1 = Recycle1 %>% left_join(anno, by = c('Run'))
  
  # Process double heavy and double Heavy-Light double lysine peptides
  recycle_LH <- as.data.frame(Recycle1 %>% filter(numb == 1))
  recycle_HH <- as.data.frame(Recycle1 %>% filter(numb == 2))
  
  sect <- intersect(recycle_LH$seqcharge,recycle_HH$seqcharge)
  recycle_LH <- recycle_LH %>% filter(seqcharge %in% sect)
  recycle_LH <- recycle_LH %>% group_by(Run,Tissue,seqcharge,Protein.Group) %>% summarise(Ms1.Area = mean(Ms1.Area,na.rm = T))
  recycle_HH <- recycle_HH %>% filter(seqcharge %in% sect)
  
  recycle_HH <- reshape2::dcast(recycle_HH, seqcharge ~ Run, value.var = 'Ms1.Area')
  rownames(recycle_HH) <- recycle_HH$seqcharge
  recycle_HH$seqcharge <- NULL
  recycle_HH <- as.matrix(recycle_HH)
  recycle_HH[log10(recycle_HH) < 6.5] <- NA
  
  recycle_LH <- reshape2::dcast(recycle_LH, seqcharge ~ Run, value.var = 'Ms1.Area')
  rownames(recycle_LH) <- recycle_LH$seqcharge
  recycle_LH$seqcharge <- NULL
  recycle_LH <- as.matrix(recycle_LH)
  recycle_LH[log10(recycle_LH) < 6.5] <- NA
  hist(log10(as.matrix(recycle_HH)))
  
  # Final estimate for frac Heavy free AA in media
  recycle_estimate_mat <- 2/(2+recycle_LH/recycle_HH)
  
  # Sample averages --  may vary for proteins as sample is bulk and different cell type specific protein expression
  # can lead to different recycling values for different proteins
  
  
  
  
  recycle_df <- as.data.frame(colMedians(as.matrix(recycle_estimate_mat),na.rm=T))
  colnames(recycle_df) <- 'Ratio'
  recycle_df$Run <- rownames(recycle_df)
  
  return(recycle_df)
  
}

# This fits the data to \gamma (t) to find the function needed downstream to find Kr
fit_sum_exp_model <- function(recycle_df){
  
  
  # IC for guess
  a_init <- log(2)/0.7
  b_init <- log(2)/20
  c_init <- .6
  
  
  # Save parameters
  age_save <- c()
  tissue_save <- c()
  a_save <- c()
  b_save <- c()
  c_save <- c()
  
  recycle_df$Age[recycle_df$Tissue == 'BM'] <- 'Old'
  
  for(i in unique(recycle_df$Age)){
    for (j in unique(recycle_df$Tissue)){
      
      recycle_df_temp <- recycle_df %>% filter(Age == i)
      recycle_df_temp <- recycle_df_temp %>% filter(Tissue == j)
      
      if(nrow(recycle_df_temp) >0){
        
        time_obs <- c(0, recycle_df_temp$Time)
        recycle_value_obs <- c(0,recycle_df_temp$Ratio)
        
        
        
        # Use nls() to fit the model and tune a and b
        fit <- nls(recycle_value_obs ~ 1 - c * exp(-a * time_obs) - c * exp(-b * time_obs), 
                   start = list(a = a_init, b = b_init,c = c_init),
                   control = nls.control(maxiter = 1000),
                   algorithm = "port")
        
        # Extract the fitted parameters
        params <- coef(fit)
        a_fitted <- params["a"]
        b_fitted <- params["b"]
        c_fitted <- params["c"]
        
        age_save <- c(age_save,i)
        tissue_save <- c(tissue_save,j)
        a_save <- c(a_save,a_fitted)
        b_save <- c(b_save,b_fitted)
        c_save <- c(c_save,c_fitted)
        
      }
    }
  }
  
  
  recycle_df_final <- as.data.frame(cbind(age_save,tissue_save,a_save,b_save,c_save))
  recycle_df_final$a_save <- as.numeric(recycle_df_final$a_save)
  recycle_df_final$b_save <- as.numeric(recycle_df_final$b_save)
  recycle_df_final$c_save <- as.numeric(recycle_df_final$c_save)
  
  
  return(recycle_df_final)
}

# Function that solves for Kr (alpha is variable used) math is in methods or can read from function
L_function <- function(alpha, L0, a, b, t) {
  c = .5
  beta = L0*alpha
  #Lt <- (exp(-(a*t) - b*t)*(-(a*beta*c*exp(a*t)) - b*beta*c*exp(b*t) + 
  #                            a*beta*c*exp(a*t + b*t) + b*beta*c*exp(a*t + b*t) + a*b*exp(a*t + b*t)*L0))/(a*b)
  
  
  Lt <- (a*beta*c - 2*alpha*beta*c + b*beta*c + alpha*beta*c*exp((-a + alpha)*t) - b*beta*c*exp((-a + alpha)*t) - 
      a*beta*c*exp((alpha - b)*t) + alpha*beta*c*exp((alpha - b)*t) - a*alpha*L0 + alpha^2*L0 + a*b*L0 - alpha*b*L0)/
    ((-a + alpha)*(alpha - b)*exp(alpha*t)) 
  
  return(Lt)
}

# RSS objective function for optim to solve for Kr
objective_function <- function(alpha, L_measured, L0, a, b, t) {
  L_predicted <- L_function(alpha, L0, a, b, t)
  residuals <- L_predicted - (L_measured)
  return(sum(residuals^2))
} 

# Run the optimization on each data point within a sample to get Kr values for a given sample
Recycle_adj <- function(L_mat, L0_mat,anno,recycle_params){
  
  
  
  alpha_recycle_adj_mat <- matrix(data= NA,nrow = nrow(L_mat),ncol = ncol(L_mat))
  
  colnames(alpha_recycle_adj_mat) <- colnames(L_mat)
  rownames(alpha_recycle_adj_mat) <- rownames(L_mat)
  
  anno$Age[anno$Tissue == 'BM'] <- 'Old'
  
  for (i in 1:ncol(L_mat)) {
    print(i)
    anno_hold <- anno %>% filter(Run == colnames(L_mat)[i])
    recycle_params_hold <- recycle_params %>% filter(tissue_save == anno_hold$Tissue & age_save == anno_hold$Age)
    t = anno_hold$Time
    
    for (j in 1:nrow(L_mat)) {
      
      # Get the current L and L0 values
      L_measured <- (L_mat[j, i])
      L0_value <- L0_mat[j, i]
      if((is.na(L0_value) + is.na(L_measured)) == 0){
        # Perform 1D optimization using "Brent" method
        result <- optim(par = .3, fn = objective_function, 
                        L_measured = L_measured, L0 = L0_value, a = recycle_params_hold$a_save, b = recycle_params_hold$b_save, t = t, 
                        method = "L-BFGS-B", lower = 0, upper = 10) # Adjust the range for alpha as needed
        
        # Store the optimal alpha in the result matrix
        alpha_recycle_adj_mat[j, i] <- result$par
      }
    }
  }
  
  return(alpha_recycle_adj_mat)

} # implement recycling


# Getting rid of impossible ratios cause presumably by noise
# by impossible we mean there isnt enough availible heavy lysine to support that ratio of heavy to light
NA_false_ratios <- function(mat,recycle_df_params,anno){
  #mat = Bulk_all_turn_pep_HL
  #recycle_df_params = recycle_df_final
  #anno = anno
  
  mat[mat==Inf] <- NA
  mat[mat==-Inf] <- NA
  mat[mat==0] <- NA
  
  
  #recycle_df$ID <- paste0(recycle_df$Time,recycle_df$Tissue,recycle_df$Age)
  anno$Age[anno$Tissue == 'BM'] <- 'Old'
  for(i in 1:ncol(mat)){
    
    Time <- anno$Time[anno$Run == colnames(mat)[i]]  
    Tissue <- anno$Tissue[anno$Run == colnames(mat)[i]]
    Age <- anno$Age[anno$Run == colnames(mat)[i]]
    
    recycle_df_temp <- recycle_df_params %>% filter(tissue_save == Tissue & age_save == Age)
    
    val <- recycle_df_temp$c_save * exp(-recycle_df_temp$a_save * Time) + recycle_df_temp$c_save * exp(-recycle_df_temp$b_save * Time)
    
    vect <- mat[,i]
    vect[vect < val] <- NA
    mat[,i] <- vect
    
  }
  
  mat <- as.matrix(mat)
  mat[mat==0]<- NA
  mat[mat==Inf]<- NA
  mat[mat==-Inf]<- NA
  
  return(mat)
} 



### Processing for bulk data

# Top 3 median is used for absolute protein quant
FindTop3 <- function(mat,PGs){
  
  for(i in 1:ncol(mat)){
    
    mat[,i] <- mat[,i]/median(mat[,i],na.rm = T)
    
  }
  
  Find_top3 <- rowMeans(mat,na.rm = T)
  Find_top3 = as.data.frame(Find_top3)
  Find_top3$prot <- PGs[,2]
  Find_top3 <- Find_top3 %>% filter(prot %in% PGsL[,2])
  Find_top3$pep <- rownames(Find_top3)
  for(i in unique(Find_top3$prot)){
    
    Find_top3_filt <- Find_top3 %>% filter(prot == i)
    
    if(nrow(Find_top3_filt) > 3){
      Find_top3_filt <- Find_top3_filt[order(-Find_top3_filt$Find_top3),]
      Find_top3_filt <- Find_top3_filt[4:nrow(Find_top3_filt),]
      Find_top3 <- Find_top3 %>% filter(!pep %in% rownames(Find_top3_filt))
    }
  }
  
  return(Find_top3)
  
} 

# collapsing peptide to protein by median top 3 most abundant
Protein_collapse <- function(mat,PGs){

  
  rownames(PGs) <-PGs[,1]
  PGs <- PGs[rownames(mat),]
  
  mat <- as.data.frame(mat)
  mat$prot <- PGs[,2]
  
  mat <- reshape2::melt(mat,ids = 'prot')
  mat <- mat %>% group_by(variable,prot) %>% summarise(value = median(value,na.rm=T))
  
  mat <- dcast(mat,prot ~ variable,value.var = 'value')
  rownames(mat) <- mat$prot
  
  
  mat$prot <- NULL
  
  mat <- as.matrix(mat)
  
  return(mat)
  
} 

# collapsing peptide to protein by median top 3 most abundant
Bulk_preproc <- function(Bulk){
  # make unique precursors
  Bulk$seqcharge <- paste0(Bulk$Stripped.Sequence,Bulk$Precursor.Charge)
  
  
  # Pulse was only with Lysine, so separate R and K peptides
  Bulk$type = str_sub(Bulk$Stripped.Sequence, -1)
  Bulk <- Bulk %>% filter(type == 'K')
  
  
  # Translation filtering for Heavy/Light peptides
  #Bulk_LysC <- Bulk_LysC %>% filter(Translated.Q.Value < .01)
  
  
  # Filter for Heavy and light
  Bulk$lab = str_sub(Bulk$Precursor.Id,-3, -3)
  Bulk_L <- Bulk %>% filter(lab == 'L')
  Bulk_H <- Bulk %>% filter(lab == 'H')
  
  
  
  # Make Heavy and Light peptide X sample matricies 
  Bulk_H <- dcast(Bulk_H,seqcharge + Protein.Group ~ Run, value.var = 'Ms1.Area')
  Bulk_L <- dcast(Bulk_L,seqcharge + Protein.Group ~ Run, value.var = 'Ms1.Area')
  
  sect = intersect(Bulk_H$seqcharge,Bulk_L$seqcharge)
  
  Bulk_H <- Bulk_H %>% filter(seqcharge %in% sect)
  Bulk_H <- Bulk_H %>% distinct(seqcharge,.keep_all = T)
  Bulk_L <- Bulk_L %>% filter(seqcharge %in% sect)
  Bulk_L <- Bulk_L %>% distinct(seqcharge,.keep_all = T)
  
  rownames(Bulk_L) <- Bulk_L$seqcharge
  rownames(Bulk_H) <- Bulk_H$seqcharge
  
  # save peptide to protein mapping
  PGsL <- cbind(Bulk_L$seqcharge,Bulk_L$Protein.Group)
  
  
  Bulk_L$Protein.Group <- NULL
  Bulk_L$seqcharge <- NULL
  Bulk_H$seqcharge <- NULL
  Bulk_H$Protein.Group <- NULL
  
  
  return(list(light = Bulk_L,heavy = Bulk_H,PGs = PGsL))
  
}

# This runs the Recycle_adj function over each sample to get (Kr X tissue sample) matrix
# !!!!! For some reason I made it log(2)/Kr.. I dont know why too lazy to change it now
# but just FYI, converting is easy anyways
Bulk_halflife_compute <- function(preproc,anno,recycle_df_final){
  Bulk_all_Kpep<- as.matrix(preproc$heavy) + as.matrix(preproc$light)
  rownames(Bulk_all_Kpep) <- preproc$PGs[,1]
  
  
  
  # Take top 3 most abundant peptides for each protein to reduce noise
  Top3_abundant <- FindTop3(Bulk_all_Kpep,preproc$PGs)
  
  
  # Data needed for model
  
  L0_mat = Bulk_all_Kpep[Top3_abundant$pep,]
  L0_mat[L0_mat==0]<- NA
  nrow(L0_mat)
  
  
  bulk_total <- L0_mat
  bulk_total_prot <- Protein_collapse(bulk_total,preproc$PGs)
  
  
  L_mat = as.matrix(preproc$light)
  rownames(L_mat) <- preproc$PGs[,1]
  L_mat[L_mat==0]<- NA
  L_mat <- L_mat[Top3_abundant$pep,]
  
  
  # Save fract Light over total and filtering out peptides with more heavy than
  # should be possible with availible AA fraction
  Bulk_all_turn_pep_HL = L_mat/L0_mat
  Bulk_all_turn_pep_HL <- NA_false_ratios(Bulk_all_turn_pep_HL,anno = anno,recycle_df_params = recycle_df_final)
  #Bulk_all_turn_pep <- -log(Bulk_all_turn_pep_HL)
  
  
  
  L0_mat[is.na(Bulk_all_turn_pep_HL)] <- NA
  L_mat[is.na(Bulk_all_turn_pep_HL)] <- NA
  
  recycle_adjusted_alpha <- Recycle_adj(L_mat,L0_mat,anno, recycle_df_final)
  halflife_pep <- log(2)/recycle_adjusted_alpha
  halflife_pep[halflife_pep > 250] <- NA
  
  halflife_prot <- Protein_collapse(halflife_pep,preproc$PGs)
  
  
  
  
  return(list(HL = halflife_prot,prot = bulk_total_prot))
}

# data processing steps ...
Bulk_proc_trypsin_abundance <- function(Bulk,Bulk_k_preproc){
  
  Bulk$seqcharge <- paste0(Bulk$Stripped.Sequence,Bulk$Precursor.Charge)
  Bulk$type = str_sub(Bulk$Stripped.Sequence, -1)
  Bulk_R <- Bulk %>% filter(type != 'K')
  
  Bulk_R$seqrun <- paste0(Bulk_R$seqcharge,Bulk_R$Run)
  Bulk_R <- Bulk_R %>% distinct(seqrun,.keep_all = T)
  Bulk_R <- reshape2::dcast(Bulk_R,seqcharge + Protein.Group ~ Run, value.var = 'Ms1.Area')
  
  PGsR <- cbind(Bulk_R$seqcharge,Bulk_R$Protein.Group)
  
  rownames(Bulk_R) <- Bulk_R$seqcharge
  Bulk_R$seqcharge <- NULL
  Bulk_R$Protein.Group <- NULL
  
  Bulk_all_Kpep <- Bulk_k_preproc$light +  Bulk_k_preproc$heavy
  
  Bulk_all_pep <- rbind(Bulk_R,Bulk_all_Kpep)
  Bulk_all_pep <- as.matrix(Bulk_all_pep)
  
  PGs <- rbind(PGsR,Bulk_k_preproc$PGs)
  rownames(Bulk_all_pep) <- PGs[,1]
  
  
  Top3_total_trypsin <- FindTop3(Bulk_all_pep,PGs)
  
  
  bulk_total_trypsin <- Bulk_all_pep[Top3_total_trypsin$pep,]
  bulk_total_trypsin[bulk_total_trypsin==Inf] <- NA
  bulk_total_trypsin[bulk_total_trypsin==-Inf] <- NA
  bulk_total_trypsin[bulk_total_trypsin==0] <- NA
  
  
  
  
  bulk_total_prot_trypsin <- Protein_collapse(bulk_total_trypsin,PGs)
  

  
  return(bulk_total_prot_trypsin)
  
}



### Plotting functions!!

plot_precursors <- function(Bulk_trypsin,Bulk_LysC){
  
  
  Bulk_trypsin$seqcharge <- paste0(Bulk_trypsin$Stripped.Sequence,Bulk_trypsin$Precursor.Charge)
  Bulk_LysC$seqcharge <- paste0(Bulk_LysC$Stripped.Sequence,Bulk_LysC$Precursor.Charge)
  
  Bulk_LysC_ <- dcast(Bulk_LysC,seqcharge~Run,value.var = 'Ms1.Area')
  Bulk_LysC_$seqcharge <- NULL
  df_pnumb <- as.data.frame(colSums(Bulk_LysC_ !=0))
  df_pnumb$Run <- rownames(df_pnumb)
  df_pnumb <- df_pnumb %>% left_join(anno, by = 'Run')
  colnames(df_pnumb)[1] <- 'Protein_numb'
  df_pnumb$dig <- 'LysC'
  
  Bulk_trypsin_ <- dcast(Bulk_trypsin,seqcharge~Run,value.var = 'Ms1.Area')
  Bulk_trypsin_$seqcharge <- NULL
  df_pnumb2 <- as.data.frame(colSums(Bulk_trypsin_ !=0))
  df_pnumb2$Run <- rownames(df_pnumb2)
  df_pnumb2 <- df_pnumb2 %>% left_join(anno, by = 'Run')
  colnames(df_pnumb2)[1] <- 'Protein_numb'
  df_pnumb2$dig <- 'Trypsin'
  
  df_pnumb <- rbind(df_pnumb2,df_pnumb)
  
  ggplot(df_pnumb, aes(x = Tissue, y = Protein_numb, color = dig)) + 
    geom_boxplot(position = position_dodge(width = 0.8), outlier.shape = NA) +  # Dodge boxplots
    geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8), 
                size = 2.5, alpha = 0.7) +  # Dodge and jitter points
    ylim(c(0, 115000)) + # Free x-scales if needed
    theme_classic(base_size = 16) +  # Set larger base font size
    theme(
      plot.title = element_text(hjust = 0.5, size = 18),  # Centered and bold title
      axis.text.x = element_text(angle = 45, hjust = 1, size = 14, color = "black"),  # Rotated x-axis labels
      axis.text.y = element_text(size = 14, color = "black"),  # Larger y-axis text
      axis.title.y = element_text(size = 16),  # Larger and bold y-axis title
      legend.title = element_text(size = 14, face = "bold"),  # Larger and bold legend title
      legend.text = element_text(size = 12),  # Larger legend text
      strip.text = element_text(size = 15, face = "bold")  # Larger facet labels
    ) +  # Use a colorblind-friendly palette
    labs(
      x = "",  # Remove x-axis title
      y = "# precursors quantified",  # Improved y-axis label
    )
}


Go_Terms_boxplot <- function(go_list,df){

  df_full_dist <- melt(df,id.vars = 'gene')
  df_full_dist$term <- 'All proteins'
  
  
  
  count <- 0
  for(i in go_list){
    
    if(count == 0){
      
      go_glu <- Go %>% filter(term == i)
      df_plt_glu <- df %>% filter(gene %in% toupper(go_glu$gene))
      df_plt_glu <- melt(df_plt_glu,id.vars = 'gene')
      df_plt_glu$term <- i
      df_add <- df_plt_glu
    }
    
    go_glu <- Go %>% filter(term == i)
    df_plt_glu <- df %>% filter(gene %in% toupper(go_glu$gene))
    df_plt_glu <- melt(df_plt_glu,id.vars = 'gene')
    df_plt_glu$term <- i
    
    df_add <- rbind(df_add,df_plt_glu)
    
    count <- 1
  }
  
  # gos <- ggplot(df_add, aes(x = term,y = as.numeric(as.character(value)), color = variable))+
  #   ylab('Fold change') +
  #   theme_classic(base_size =18) + xlab('') + ylim(c(-3,3))+
  #   geom_boxplot(outlier.shape=NA, width=0.4, position = position_dodge(width=0.6)) +
  #   geom_point(aes(color = variable), position=position_jitterdodge(dodge.width = 0.6, jitter.width = 0.1), alpha=1) 
  # 
  # full <- ggplot(df_full_dist, aes(x = term,y = as.numeric(as.character(value)), color = variable))+
  #   ylab('Fold change') +
  #   theme_classic(base_size =18) + xlab('') + ylim(c(-3,3))+
  #   geom_boxplot(width=0.4, position = position_dodge(width=0.6))
  # 
  # full+gos
  
  combined_plot <- ggplot() +
    # Boxplot for the full dataset
    geom_boxplot(
      data = df_full_dist,
      aes(x = term, y = as.numeric(as.character(value)), color = variable),
      width = 0.4,
      position = position_dodge(width = 0.6),
      outlier.shape = NA  # Suppress outliers in the full dataset
    ) +
    # Boxplot for the subset
    geom_boxplot(
      data = df_add,
      aes(x = term, y = as.numeric(as.character(value)), color = variable),
      width = 0.3,
      position = position_dodge(width = 0.6),
      outlier.shape = NA  # Suppress outliers in the subset
    ) +
    # Jittered points for the subset
    geom_point(
      data = df_add,
      aes(x = term, y = as.numeric(as.character(value)), color = variable),
      position = position_jitterdodge(dodge.width = 0.6, jitter.width = 0.1),
      size = 1.5,
      alpha = 0.8  # Adjust transparency
    ) +
    facet_wrap(~term, scales = "free_x", nrow = 1) +  # Separate by term
    ylab('Fold change') +
    xlab('') +
    ylim(c(-3, 3)) +
    theme_classic(base_size = 16) +
    geom_hline(yintercept = 0) +
    theme(axis.text.x = element_text(angle = 45,hjust=1,vjust=1))+
    scale_color_manual(values = c('#F3766E','#000000'))
  
  # Display the plot
  combined_plot
  
}


plot_prot_cor <- function(prot_list,anno,df_abs1,df_hl1,df_abs2,df_hl2){
  
  
  for(i in 1:length(prot_list)){
    if(i==1){
      df1 <- as.data.frame(cbind(df_hl1[prot_list[i],],df_abs1[prot_list[i],]))
      colnames(df1) <- c('HL','Abundance')
      df1$prot <- prot_list[i]
      df1$Sample <- colnames(df_abs1)
      df1$digest = 'Trypsin'
      df2 <- as.data.frame(cbind(df_hl2[prot_list[i],],df_abs2[prot_list[i],]))
      colnames(df2) <- c('HL','Abundance')
      df2$prot <- prot_list[i]
      df2$Sample <- colnames(df_abs2)
      df2$digest = 'LysC'
      
      df <- rbind(df1,df2)
      
      
    }else{
      df1 <- as.data.frame(cbind(df_hl1[prot_list[i],],df_abs1[prot_list[i],]))
      colnames(df1) <- c('HL','Abundance')
      df1$prot <- prot_list[i]
      df1$Sample <- colnames(df_abs1)
      df1$digest = 'Trypsin'
      df2 <- as.data.frame(cbind(df_hl2[prot_list[i],],df_abs2[prot_list[i],]))
      colnames(df2) <- c('HL','Abundance')
      df2$prot <- prot_list[i]
      df2$Sample <- colnames(df_abs2)
      df2$digest = 'LysC'
      df1 <- rbind(df1,df2)
      
      df <- rbind(df,df1)
      
    }
    

    
  }

  r2_values <- df %>%
    group_by(prot) %>%
    summarize(
      r2 = round(cor(HL, Abundance, use = "complete.obs")^2, 2)  # Compute R^2
    )
  
  df <- df %>% left_join(anno, by = 'Sample')
  
  # Plot
  ggplot(df, aes(x = log2(log(2)/2^HL), y = Abundance,shape = digest,color = Tissue)) +
    geom_point() +
    #xlim(c(-3, 2)) +
    #ylim(c(-3, 2)) +
    facet_wrap(~prot,nrow = 1) +
    geom_text(
      data = r2_values,  # Use summarized data for labels
      aes(label = paste0("R^2 = ", r2), x = 2, y = 2),  # Position for the text
      inherit.aes = FALSE,  # Prevent inheriting aesthetics from ggplot
      hjust = 1, vjust = 1.5
    ) +
    theme_bw(base_size = 20) +
    labs(x = "Half-Life (HL)", y = "Abundance") +
    theme(strip.text.x = element_text(size = 14))+
    scale_color_manual(values = c('#BC3EFA','#088109','#FB8A02','#0237FB'))
  
}


Go_Terms_boxplot_cor <- function(go_list,df){
  df$prot <- NULL
  df_full_dist <- df
  df_full_dist$term <- 'All proteins'
  df_full_dist$split_gene <- NULL
  df_full_dist <- melt(df_full_dist,id.vars = 'term')
  
  df$null_dist <- NULL
  
  count <- 0
  for(i in go_list){
    
    if(count == 0){
      
      go_glu <- Go %>% filter(term == i)
      df_plt_glu <- df %>% filter(split_gene %in% go_glu$gene)
      df_plt_glu$term <- i
      df_add <- df_plt_glu
    }
    
    go_glu <- Go %>% filter(term == i)
    df_plt_glu <- df %>% filter(split_gene %in% go_glu$gene)
    df_plt_glu$term <- i
    df_add <- rbind(df_add,df_plt_glu)
    count <- 1
  }

  
  combined_plot <- ggplot() +
    # Boxplot for the full dataset
    geom_boxplot(
      data = df_full_dist,
      aes(x = variable, y = as.numeric(as.character(value))),
      width = .5,
      outlier.shape = NA
        # Suppress outliers in the full dataset
    ) +
    # Boxplot for the subset
    geom_boxplot(
      data = df_add,
      aes(x = term, y = cor),
      width = .5,
      outlier.shape = NA
        # Suppress outliers in the subset
    ) +
    # Jittered points for the subset
    geom_jitter(
      data = df_add,
      aes(x = term, y = as.numeric(as.character(cor))),
      width = .2,
      alpha = 0.5  # Adjust transparency
    ) +
    facet_wrap(~term, scales = "free_x", nrow = 1) +  # Separate by term
    ylab('Fold change') +
    xlab('') +
    ylim(c(-.5, 1)) +
    theme_classic(base_size = 16) +
    geom_hline(yintercept = 0) +
    theme(axis.text.x = element_text(angle = 45,hjust=1,vjust=1))
  
  # Display the plot
  combined_plot
  
}

