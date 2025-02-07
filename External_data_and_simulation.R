# You have to change the path name here 
source('~/Desktop/Projects/Miceotopes/Bulk/Functions_for_analysis.R')



###### Supp Figure 1 log normal

#Schwan <- read.csv('/Users/andrewleduc/Desktop/Projects/Dilution_paper/External/Datasets/Schwanhauser_et_al.csv')
Schwan <- read.csv('https://zenodo.org/records/14827610/files/Schwanhauser_et_al.csv?download=1')
colnames(Schwan) <- c('Protein','Abundance','half.life')


# Abundance
set.seed(123)
X <- sample(log10(Schwan$Abundance),3000)  # Replace this with your actual vector

# Perform a Shapiro-Wilk Normality Test
shapiro_test <- shapiro.test(X)
cat("Shapiro-Wilk Test:\n")
print(shapiro_test)

# Create a data frame for ggplot
data <- data.frame(X = X)

# Generate a QQ plot with ggplot2
ggplot(data, aes(sample = X)) +
  stat_qq(size = 2, color = "blue") +                  # Add QQ plot points
  stat_qq_line(color = "red", linetype = "dashed") +   # Add reference line
  labs(
    title = "log Protein abundance, p-value = 1.03e-5",
    x = "Theoretical Quantiles",
    y = "Sample Quantiles"
  ) +
  theme_minimal(base_size = 14)


mean_X <- mean(X)
sd_X <- sd(X)

# Create the density plot with the Gaussian approximation
ggplot(data, aes(x = X)) +
  geom_histogram(aes(y = ..density..), bins = 30, fill = "skyblue", color = "black", alpha = 0.7) +  # Histogram
  stat_function(fun = dnorm, args = list(mean = mean_X, sd = sd_X), 
                color = "red", size = 1, linetype = "dashed") +  # Gaussian curve
  labs(
    title = "Distribution of Data with Gaussian Approximation",
    x = "Log 10 copy number",
    y = "# Proteins"
  ) +
  theme_minimal(base_size = 14)


# Deg rates

set.seed(123)
X <- sample(log2(Schwan$half.life),3000)  # Replace this with your actual vector

# Perform a Shapiro-Wilk Normality Test
shapiro_test <- shapiro.test(X)
cat("Shapiro-Wilk Test:\n")
print(shapiro_test)

ks_test <- ks.test(X, "pnorm", mean = mean(X), sd = sd(X))
cat("Kolmogorov-Smirnov Test:\n")
print(ks_test)

# Create a data frame for ggplot
data <- data.frame(X = X)

# Generate a QQ plot with ggplot2
ggplot(data, aes(sample = X)) +
  stat_qq(size = 2, color = "blue") +                  # Add QQ plot points
  stat_qq_line(color = "red", linetype = "dashed") +   # Add reference line
  labs(
    title = "log Half-life, p-value = 8.26e-10",
    x = "Theoretical Quantiles",
    y = "Sample Quantiles"
  ) +
  theme_minimal(base_size = 14)

X <- sample(log2(Schwan$half.life))  # Replace this with your actual vector
data <- data.frame(X = X)
X = log2(log(2)/2^X)
data$X <- log2(log(2)/2^data$X)

mean_X <- mean(X)
sd_X <- sd(X)

# Create the density plot with the Gaussian approximation
ggplot(data, aes(x = X)) +
  geom_histogram(aes(y = ..density..), bins = 40, fill = "skyblue", color = "black", alpha = 0.7) +  # Histogram
  stat_function(fun = dnorm, args = list(mean = mean_X, sd = sd_X), 
                color = "red", size = 1, linetype = "dashed") +  # Gaussian curve
  labs(
    title = "Distribution of Data with Gaussian Approximation",
    x = "Log 2 half life (hours)",
    y = "# Proteins"
  ) +
  theme_minimal(base_size = 14)+
  ylim(c(0,.4))



# Translation rates

set.seed(123)
X <- sample(log2(Schwan$translation),3000)  # Replace this with your actual vector

# Perform a Shapiro-Wilk Normality Test
shapiro_test <- shapiro.test(X)
cat("Shapiro-Wilk Test:\n")
print(shapiro_test)

# Create a data frame for ggplot
data <- data.frame(X = X)

# Generate a QQ plot with ggplot2
ggplot(data, aes(sample = X)) +
  stat_qq(size = 2, color = "blue") +                  # Add QQ plot points
  stat_qq_line(color = "red", linetype = "dashed") +   # Add reference line
  labs(
    title = "log Translation rate, p-value = 8.26e-10",
    x = "Theoretical Quantiles",
    y = "Sample Quantiles"
  ) +
  theme_minimal(base_size = 14)

mean_X <- mean(X)
sd_X <- sd(X)

# Create the density plot with the Gaussian approximation
ggplot(data, aes(x = X)) +
  geom_histogram(aes(y = ..density..), bins = 30, fill = "skyblue", color = "black", alpha = 0.7) +  # Histogram
  stat_function(fun = dnorm, args = list(mean = mean_X, sd = sd_X), 
                color = "red", size = 1, linetype = "dashed") +  # Gaussian curve
  labs(
    title = "Distribution of Data with Gaussian Approximation",
    x = "Log 10 Translation rate (copy/hr)",
    y = "# Proteins"
  ) +
  theme_minimal(base_size = 14)



###### Supp Figure 1 Translation vs Degradation co-linearity

#Sav_Bcell <- read.csv('/Users/andrewleduc/Desktop/Projects/Dilution_paper/External/PythonProcData/Bcells.csv')
Sav_Bcell <- read.csv('https://zenodo.org/records/14827610/files/Bcells.csv?download=1')
Sav_Bcell$Trans1 <- Sav_Bcell$Abundance1 + Sav_Bcell$Deg_rate1
b_cor <- cor(Sav_Bcell$Trans1,Sav_Bcell$Deg_rate1)
#plot(Sav_Bcell$Trans1,Sav_Bcell$Deg_rate1)


#Sav_Hep <- read.csv('/Users/andrewleduc/Desktop/Projects/Dilution_paper/External/PythonProcData/Hepatocytes.csv')
Sav_Hep <- read.csv('https://zenodo.org/records/14827610/files/Hepatocytes.csv?download=1')
Sav_Hep$Trans1 <- Sav_Hep$Abundance1 + Sav_Hep$Deg_rate1
H_cor <- cor(Sav_Hep$Trans1,Sav_Hep$Deg_rate1)
#plot(Sav_Hep$Trans1,Sav_Hep$Deg_rate1)

#Sav_Mon <- read.csv('/Users/andrewleduc/Desktop/Projects/Dilution_paper/External/PythonProcData/Monocytes.csv')
Sav_Mon <- read.csv('https://zenodo.org/records/14827610/files/Monocytes.csv?download=1')
Sav_Mon$Trans1 <- Sav_Mon$Abundance1 + Sav_Mon$Deg_rate1
Mon_cor <- cor(Sav_Mon$Trans1,Sav_Mon$Deg_rate1)
#plot(Sav_Mon$Trans1,Sav_Mon$Deg_rate1)

#Sav_NK <- read.csv('/Users/andrewleduc/Desktop/Projects/Dilution_paper/External/PythonProcData/NKcells.csv')
Sav_NK <- read.csv('https://zenodo.org/records/14827610/files/NKcells.csv?download=1')
Sav_NK$Trans1 <- Sav_NK$Abundance1 + Sav_NK$Deg_rate1
NK_cor <- cor(Sav_NK$Trans1,Sav_NK$Deg_rate1)
#plot(Sav_NK$Trans1,Sav_NK$Deg_rate1)

Schwan$translation <- Schwan$Abundance * (log(2)/Schwan$half.life + log(2)/24)




cor_store <- c(NK_cor,Mon_cor,H_cor,b_cor,cor(log2(Schwan$translation),log2(log(2)/Schwan$half.life)))
sample_type <- c('NK cells','Monocytes','Hepatocytes','B Cells','Mouse embryonic fibroblast')
data_set <- c('No growth','No growth','No growth','No growth','Growing')
df_plot <- as.data.frame(cbind(cor_store,sample_type,data_set))
df_plot$cor_store <- as.numeric(df_plot$cor_store )

ggplot(df_plot, aes(x = sample_type,y = cor_store)) + geom_point(size = 5) + 
  facet_wrap(~data_set,scales = 'free_x') +
  theme_bw(base_size = 15) +  # Classic theme with larger font
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),  # Centered, bold title
    axis.text = element_text(color = "black"),  # Black axis text
    axis.title = element_text(face = "bold")  # Bold axis labels
  )+
  xlab('') + ylab('Cor(log Translation x log Degradation) ') + 
  ylim(c(-.5,.5))



Sav_NK2 <- Sav_NK %>% group_by(Protein) %>% summarise(Abundance1 = median(Abundance1,na.rm=T),
                                                      Abundance2 = median(Abundance2,na.rm=T),
                                                      Deg_rate1 = median(Deg_rate1,na.rm=T),
                                                      Deg_rate2 = median(Deg_rate2,na.rm=T))

Sav_NK2$trans <- Sav_NK2$Abundance1+Sav_NK2$Deg_rate1
sd(Sav_NK2$Deg_rate1)
sd(Sav_NK2$trans)

compute_rsq(Sav_NK2$Abundance1,Sav_NK2$Deg_rate1,Sav_NK2$trans)

Sav_Mon2 <- Sav_Mon %>% group_by(Protein) %>% summarise(Abundance1 = median(Abundance1,na.rm=T),
                                                        Abundance2 = median(Abundance2,na.rm=T),
                                                      Deg_rate1 = median(Deg_rate1,na.rm=T),
                                                      Deg_rate2 = median(Deg_rate2,na.rm=T))

Sav_Mon2$trans <- Sav_Mon2$Abundance1+Sav_Mon2$Deg_rate1
sd(Sav_Mon2$Deg_rate1)
sd(Sav_Mon2$trans)


compute_rsq(Sav_Mon2$Abundance1,Sav_Mon2$Deg_rate1,Sav_Mon2$trans)

Sav_Bcell2 <- Sav_Bcell %>% group_by(Protein) %>% summarise(Abundance1 = median(Abundance1,na.rm=T),
                                                        Abundance2 = median(Abundance2,na.rm=T),
                                                        Deg_rate1 = median(Deg_rate1,na.rm=T),
                                                        Deg_rate2 = median(Deg_rate2,na.rm=T))

Sav_Bcell2$trans <- Sav_Bcell2$Abundance1+Sav_Bcell2$Deg_rate1
sd(Sav_Bcell2$Deg_rate1)
sd(Sav_Bcell2$trans)




compute_rsq(Sav_Bcell2$Abundance1,Sav_Bcell2$Deg_rate1,Sav_Bcell2$trans)

Sav_Hep2 <- Sav_Hep %>% group_by(Protein) %>% summarise(Abundance1 = median(Abundance1,na.rm=T),
                                                        Abundance2 = median(Abundance2,na.rm=T),
                                                        Deg_rate1 = median(Deg_rate1,na.rm=T),
                                                        Deg_rate2 = median(Deg_rate2,na.rm=T))


Sav_Hep2$trans <- Sav_Hep2$Abundance1+Sav_Hep2$Deg_rate1
sd(Sav_Hep2$Deg_rate1)
sd(Sav_Hep2$trans)

compute_rsq(Sav_Hep2$Abundance1,Sav_Hep2$Deg_rate1,Sav_Hep2$trans)


sd_kd <- c(sd(Sav_Hep2$Deg_rate1),sd(Sav_Bcell2$Deg_rate1),sd(Sav_Mon2$Deg_rate1),sd(Sav_NK2$Deg_rate1),sd(log2(log(2)/Schwan$half.life)))
sd_s <- c(sd(Sav_Hep2$trans),sd(Sav_Bcell2$trans),sd(Sav_Mon2$trans),sd(Sav_NK2$trans),sd(log2(Schwan$translation)))
cts <- c('Hepatocyte','B Cells','Monocytes','NK cells','Fibroblasts')

plt_sd <- data.frame(sd_kd = sd_kd, sd_s = sd_s,cts=cts)

ggplot(plt_sd,aes(y = sd_s,x = sd_kd,color = cts)) + geom_point(size=6 ) +
  dot_plot 

###### Supp Figure 1 Cor^2 vs R^2 for increasing colinearity

R_sq_save <- c()
Cor_sq_save <- c()
colinearity <- c()
slope_save <- c()


CC = 0
for(i in 1:95){
  
  
  mu_logP <- mean(log2(Schwan$translation))
  mu_logAlpha <- mean(log2(log(2)/Schwan$half.life)) 
  sigma_logP <- sd(log2(Schwan$translation))
  sigma_logAlpha <- sd(log2(log(2)/Schwan$half.life))
  correlation <- -.24  # Given correlation
  
  # Covariance matrix
  cov_matrix <- matrix(c(sigma_logP^2, correlation * sigma_logP * sigma_logAlpha, 
                         correlation * sigma_logP * sigma_logAlpha, sigma_logAlpha^2), nrow = 2)
  
  
  simulated_data <- mvrnorm(n = 20000, mu = c(mu_logP, mu_logAlpha), Sigma = cov_matrix)
  
  
  deg <- (simulated_data[, 2])
  trans <- (simulated_data[, 1])
  

  
  P = trans - deg
  cor(P,deg)

  rsq <- compute_rsq(P,deg,trans)
  sloper <- abs(TLS(deg,P)[[1]])
  
  R_sq_save <- c(R_sq_save,rsq)
  Cor_sq_save <- c(Cor_sq_save,cor(deg,P)^2)
  slope_save <- c(slope_save,sloper)
  colinearity <- c(colinearity,CC)
  
  CC = CC - .01
  
}

df_supp1c <- as.data.frame(cbind(R_sq_save,Cor_sq_save,slope_save,colinearity))
df_supp1c <- reshape2::melt(df_supp1c,id.var = 'colinearity')

ggplot(df_supp1c, aes(x = -colinearity,y = value, color = variable)) + geom_smooth()+
  dot_plot+ xlab('Cor(k_s,k_r)') +
  ggtitle('') + ylab('Removal X Abundance')



###### Figure 1 d
Schwan$translation <- Schwan$Abundance*(log(2)/Schwan$half.life)


cor_save <- c()
sd_save <- c()
dilute_save <- c()
i_save <- c()

x_range <-  x_range <- 10^(seq(log10(1), log10(1000), length.out = 100)) #seq(1, 1000, length.out = 100) 
y_range <- seq(0, 3, length.out = 150)   # Linear scale


sd_add = 0
for(i in 1:150){
  
  mu_logP <- mean(log2(Schwan$translation))
  mu_logAlpha <- mean(log2(log(2)/Schwan$half.life)) 
  sigma_logP <- sd(log2(Schwan$translation))
  sigma_logAlpha <- sd(log2(log(2)/Schwan$half.life)) + y_range[i]
  correlation <- -.15  # Given correlation
  
  # Covariance matrix
  cov_matrix <- matrix(c(sigma_logP^2, correlation * sigma_logP * sigma_logAlpha, 
                         correlation * sigma_logP * sigma_logAlpha, sigma_logAlpha^2), nrow = 2)
  
  
  simulated_data <- mvrnorm(n = 20000, mu = c(mu_logP, mu_logAlpha), Sigma = cov_matrix)
  
  
  deg <- (simulated_data[, 2])
  trans <- (simulated_data[, 1])
  
  #double_time <- 1000
  for(j in 1:101){
    
    if(j != 101){
      double_time <- x_range[j]
    }else{
      double_time <- 100000
      print('Here')
    }
    
    deg_plus_dilution <- log2(2^deg +log(2)/double_time)
    
    
    
    P = trans - deg_plus_dilution
    
    
    rsq <- compute_rsq(P,deg_plus_dilution,trans)
    
    
    cor_save <- c(cor_save, rsq)
    dilute_save <- c(dilute_save,double_time)
    sd_save <- c(sd_save,sigma_logAlpha/sigma_logP)#sd(deg_plus_dilution)/sigma_logP)
    i_save <- c(i_save,i)
    #double_time = double_time -10
    
  }  
    
  sd_add = sd_add + .2

}

data <- data.frame(X = dilute_save,Y = sd_save, Z = cor_save)

data1 <- data %>% filter(X != 100000)

# Create the plot
ggplot(data1, aes(x = X, y = (Y), fill =  (Z))) +
  geom_raster() +  # Use geom_tile() for grid data
  scale_fill_gradient2(
    low = "white",  # Low values
    high = "black",  # Mid value (Maroon)
    mid = "red",  # High values
    midpoint = 0.25,  # Define the midpoint for the gradient
    name = "R-squared"
  )+   # Adjust gradient
  labs(
    x = "Cell Division Time (h)",
    y = "Sigma_s/Sigma_r",
    fill = "R-squared"
  ) +
  theme_classic()  +
  coord_cartesian(xlim = c(2,780),ylim = c(.52,1.4)) +
  dot_plot+ scale_x_log10()


# Create the plot
data2 <- data %>% filter(X == 100000)
data2$X <- as.character(data2$X)
ggplot(data2, aes(x = X, y = (Y), fill =  (Z))) +
  geom_raster() +  # Use geom_tile() for grid data
  scale_fill_gradient2(
    low = "white",  # Low values
    high = "black",  # Mid value (Maroon)
    mid = "red",  # High values
    midpoint = 0.25,  # Define the midpoint for the gradient
    name = ""
  )+   # Adjust gradient
  labs(
    x = "",
    y = "",
    fill = ""
  ) +
  theme_classic() + 
  coord_cartesian(ylim = c(.52,1.4)) +
  dot_plot




Schwan$translation <- log(2)/Schwan$half.life*Schwan$Abundance
Schwan$translation <- Schwan$Abundance * (log(2)/Schwan$half.life + log(2)/24)

mu_logP <- mean(log2(Schwan$translation))
mu_logAlpha <- mean(log2(log(2)/Schwan$half.life + log(2)/17)) 
sigma_logP <- sd(log2(Schwan$translation))
sigma_logAlpha <- sd(log2(log(2)/Schwan$half.life + log(2)/17)) 
correlation <- -.22  # Given correlation

# Covariance matrix
cov_matrix <- matrix(c(sigma_logP^2, correlation * sigma_logP * sigma_logAlpha, 
                       correlation * sigma_logP * sigma_logAlpha, sigma_logAlpha^2), nrow = 2)


simulated_data <- mvrnorm(n = 20000, mu = c(mu_logP, mu_logAlpha), Sigma = cov_matrix)


deg <- (simulated_data[, 2])
trans <- (simulated_data[, 1])

P = trans - deg

compute_rsq(P,deg,trans)


P <- P-mean(P)
deg <- deg-mean(deg)

take <- sample(1:length(P),500)
deg <- deg[take]
P <- P[take]

data <- data.frame(
  deg = deg,  # Example values for -deg
  P = P     # Example values for P
)

# Create the ggplot
ggplot(data, aes(x = P, y = -deg)) +
  geom_point(size = 3, color = "black") +  # Data points
  geom_abline(slope = 1, intercept = 0, linetype = "solid", color = "black") +  # Line y = x
  geom_segment(aes(x = P, y = -deg, xend = P, yend = P), linetype = "dashed", color = "gray") +  # Vertical dashed lines
  annotate("text", x = max(data$P) * 0.8, y = max(-data$deg) * 0.8, 
           label = "Rsq = 0.10", size = 5, color = "black", hjust = 0) +  # Add text annotation
  labs(
    x = "P",
    y = "-K_d",
    title = "Evaluate contribution of degradation rates to protein levels"
  ) +
  dot_plot


############ Figure 1b ##############


b_rsq <- compute_rsq(Sav_Bcell2$Abundance1,Sav_Bcell2$Deg_rate1,Sav_Bcell2$trans)
b_sd <- sd(Sav_Bcell2$Deg_rate1)/sd(Sav_Bcell2$trans)

H_rsq <- compute_rsq(Sav_Hep2$Abundance1,Sav_Hep2$Deg_rate1,Sav_Hep2$trans)
H_sd <- sd(Sav_Hep2$Deg_rate1)/sd(Sav_Hep2$trans)

Mon_rsq <- compute_rsq(Sav_Mon2$Abundance1,Sav_Mon2$Deg_rate1,Sav_Mon2$trans)
Mon_sd <- sd(Sav_Mon$Deg_rate1)/sd(Sav_Mon$trans)

NK_rsq <- compute_rsq(Sav_NK2$Abundance1,Sav_NK2$Deg_rate1,Sav_NK2$trans)
NK_sd <- sd(Sav_NK2$Deg_rate1)/sd(Sav_NK2$trans)
NK_rsq <- cor(Sav_NK2$Abundance1,Sav_NK2$Deg_rate1)^2

ggplot(Sav_NK2, aes(x = 2^Abundance1, y = 2^Deg_rate1)) +
  ggpointdensity::geom_pointdensity() +
  scale_color_viridis_c() +  # Adjust line thickness for visibility
  ylab("Protein abundance") +  # Updated x-axis label
  xlab(expression("(" * K[d] + k[g] * "), 1/hour")) +  # Use expression for better formatting of subscripts
  theme_classic(base_size = 15) +  # Increase base font size
  theme(
    axis.text = element_text(size = 14, color = "black"),  # Larger axis text
    axis.title = element_text(size = 16),  # Bold axis titles
    panel.border = element_rect(color = "black", fill = NA, size = 1.5),  # Add border of consistent thickness
    axis.ticks = element_line(size = 1.5),  # Make axis ticks the same thickness
    axis.line = element_blank(),  # Remove separate axis lines to use the border
    legend.position = "none"  # Remove the legend
  ) +
  ggtitle("NK Cells") +
  scale_y_log10(labels = trans_format("log10", math_format(10^.x))) +
  scale_x_log10(labels = trans_format("log10", math_format(10^.x)))




####################################################################################
# Figure 1c 
####################################################################################

Schwan$translation <- (log(2)/Schwan$half.life+log(2)/17)*Schwan$Abundance


compute_rsq(log2(Schwan$Abundance),log2(log(2)/Schwan$half.life+log(2)/17),log2(Schwan$translation))

cor_store <- c(NK_rsq,Mon_rsq,H_rsq,b_rsq,compute_rsq(log2(Schwan$Abundance),log2(log(2)/Schwan$half.life+log(2)/17),log2(Schwan$translation)))
sample_type <- c('NK cells',' Monocytes','Hepatocytes','B Cells','Fibroblast')
data_set <- c('No growth','No growth','No growth','No growth','Growing')
df_plot <- as.data.frame(cbind(cor_store,sample_type,data_set))
df_plot$cor_store <- as.numeric(df_plot$cor_store )
df_plot$cor_store <- df_plot$cor_store * 100
ggplot(df_plot, aes(x = sample_type,y = cor_store)) + geom_point(size = 5) + 
  facet_wrap(~data_set,scales = 'free_x') +
  theme_bw(base_size = 15) +  # Classic theme with larger font
  theme(
    plot.title = element_text(hjust = 0.5),  # Centered, bold title
    axis.text = element_text(color = "black"),  # Black axis text
    axis.title = element_text(face = "bold")  # Bold axis labels
  )+
  xlab('') + ylab('Variance explained, %') + 
  ylim(c(0,60)) 






##################################################################################################
# Figure 2 analysis
##################################################################################################
Immune_act <- read.csv("https://zenodo.org/records/14827610/files/Act_rest_Bcells.csv?download=1")
#Immune_act <- read.csv('/Users/andrewleduc/Desktop/Projects/Dilution_paper/External/Datasets/Act_rest_immune.csv')
hLs <- read.csv('https://zenodo.org/records/14827610/files/Savitski_halflives.csv?download=1')
#hLs <- read.csv('/Users/andrewleduc/Desktop/Projects/Dilution_paper/External/Datasets/Savitski_halflives.csv')

# Filter data
hLs <- hLs %>% filter(Bcells.replicate.1.half_life != "#N/A")
hLs <- hLs %>% filter(Bcells.replicate.2.half_life != "#N/A")
hLs <- hLs %>% filter(Bcells.replicate.1.half_life != "Inf" )
hLs <- hLs %>% filter(Bcells.replicate.2.half_life != "Inf" )
hLs <- hLs %>% filter(Bcells.replicate.1.half_life > .9)
hLs <- hLs %>% filter(Bcells.replicate.1.dataQual =='good')

# Intersect proteins
Immune_act <- Immune_act %>% filter(Gene.names %in% hLs$gene_name)
hLs <- hLs %>% filter(gene_name %in% Immune_act$Gene.names)
hLs <- hLs %>% distinct(gene_name,.keep_all = T)
Immune_act <- Immune_act %>% distinct(Gene.names,.keep_all = T)
rownames(hLs) <- hLs$gene_name 
hLs <- hLs[Immune_act$Gene.names,]
hLs$gene_name == Immune_act$Gene.names


# fix variables
hLs$Bcells.replicate.1.half_life <- as.numeric(hLs$Bcells.replicate.1.half_life)
hLs$Bcells.replicate.2.half_life <- as.numeric(hLs$Bcells.replicate.2.half_life)

# Take averages for all memory and naive Bcells because the savitski data uses both
vect <- log2(rowMeans(cbind(Immune_act$Intensity_B.memory_02_activated,Immune_act$Intensity_B.memory_03_activated,Immune_act$Intensity_B.naive_01_activated,Immune_act$Intensity_B.naive_02_activated),na.rm=T)/
               rowMeans(cbind(Immune_act$Intensity_B.memory_02_steady.state,Immune_act$Intensity_B.memory_03_steady.state,Immune_act$Intensity_B.naive_01_steady.state,Immune_act$Intensity_B.naive_02_steady.state),na.rm=T))
vect[abs(vect) > 10] <- NA

#Half lives average from savitsski replicates
vect2 <- rowMeans(cbind(hLs$Bcells.replicate.1.half_life,hLs$Bcells.replicate.2.half_life),na.rm=T)

cor(vect,log2(vect2),use = 'pairwise.complete.obs')^2
plot(vect,log2(vect2))


# Get the reliability by comparing replicate measurments for concentrations
v1 <- log2(rowMeans(cbind(Immune_act$Intensity_B.memory_03_activated,Immune_act$Intensity_B.naive_02_activated),na.rm=T))
v1 = v1 - median(v1,na.rm=T)
v2 <- log2(rowMeans(cbind(Immune_act$Intensity_B.memory_02_activated,Immune_act$Intensity_B.naive_01_activated),na.rm=T))
v2 = v2 - median(v2,na.rm=T)
v1[abs(v1) > 10] <- NA
v2[abs(v2) > 10] <- NA

v1 <- log2(Immune_act$Intensity_B.memory_02_activated/Immune_act$Intensity_B.memory_02_steady.state) - median(log2(Immune_act$Intensity_B.memory_02_activated/Immune_act$Intensity_B.memory_02_steady.state),na.rm = T)
v2 <- log2(Immune_act$Intensity_B.memory_03_activated/Immune_act$Intensity_B.memory_03_steady.state) - median(log2(Immune_act$Intensity_B.memory_03_activated/Immune_act$Intensity_B.memory_03_steady.state),na.rm = T)

v1[v1 == Inf] <- NA
v2[v2 == Inf] <- NA
v1[v1 == -Inf] <- NA
v2[v2 == -Inf] <- NA

cor(v1,v2,use = 'pairwise.complete.obs')
plot(v1,v2, xlab = 'Replicate 1', ylab = 'Replicate 2', main = 'Cor = 0.61')

# Compute Half life reliability
cor(log2(hLs$Bcells.replicate.1.half_life),log2(hLs$Bcells.replicate.2.half_life),use = 'pairwise.complete.obs')
plot(log2(hLs$Bcells.replicate.1.half_life),log2(hLs$Bcells.replicate.2.half_life), xlab = 'Replicate 1', ylab = 'Replicate 2', main = 'Cor = 0.95')


cor(vect,log2(vect2),use = 'pairwise.complete.obs')^2/(.95*cor(v1,v2,use = 'pairwise.complete.obs'))
plot(vect-median(vect,na.rm=T),log2(vect2),xlab = 'log2(Resting/Activated) - Abundance',ylab = 'log2(Half life) hours',
     main = 'Raw Rsq = 0.21, Noise corrected Rsq = 0.35')

##################################################################################################
# Figure 2 b
##################################################################################################

df_rank <- as.data.frame(cbind(log2(log(2)/vect2),log2(log(2)/vect2 + log(2)/16))) 
df_rank <- df_rank[order(df_rank$V1),]
df_rank$rank  <- 1:nrow(df_rank)
df_rank <- melt(df_rank,id.vars = 'rank')

ggplot(df_rank, aes(y = 2^value, x = rank, color = variable )) + geom_point()+
  xlab('Rank') + ylab('K_r')+theme_bw(base_size = 15) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),  # Centered, bold title
    axis.text.y = element_text(color = "black",size = 18)  # Black axis text
  )+
    scale_y_log10(
      labels = trans_format("log10", math_format(10^.x)), 
      breaks = scales::trans_breaks("log10", function(x) 10^x)  # Custom breaks
    ) +
  scale_color_manual(values = c('#7F3F98','#3FB655'))
    
  



##################################################################################################
# Figure 2 c
##################################################################################################

df_plt <- as.data.frame(cbind(vect,log2(log(2)/vect2) - log2(log(2)/vect2 + log(2)/16))) 
df_plt$gene <- hLs$gene_name
df_plt$col <- 'normal'
df_plt$col[df_plt$gene %in% c('HIST2H2AB','HIST2H3PS2','HIST1H4A','HIST1H1C','HIST1H1E','HIST1H1B','HIST1H1D')] <- 'Histone'
df_plt$dif = (df_plt$vect - df_plt$V2)
df_plt$HL <- vect2

cor(log2(df_plt$HL),df_plt$dif, use = 'pairwise.complete.obs')

library(scales)
ggplot(df_plt,aes(y = 2^(vect - median(vect,na.rm=T)),x = 2^(V2- median(V2,na.rm=T)),color = col,alpha = col)) + geom_point(size = 3)+ 
  xlab('Expected dilution') + ylab('Concentration change') +
  geom_abline(intercept = 0, slope = 1)+
  dot_plot + scale_x_log10(labels = trans_format("log10", math_format(10^.x))) + 
  scale_y_log10(labels = trans_format("log10", math_format(10^.x))) +
  scale_color_manual(values = c('black','grey50'))+
  scale_alpha_manual(values = c(1,.4))


##################################################################################################
# Figure 2 d
##################################################################################################


df_plt$regulated <- (df_plt$vect - median(df_plt$vect,na.rm=T)) -  (df_plt$V2- median(df_plt$V2,na.rm=T))


df_plt2 <- df_plt[order(-df_plt$regulated),]
df_plt2$rank <- 1:nrow(df_plt2)


ggplot(df_plt,aes(y = 2^regulated,x = HL,color = col,alpha = col)) + geom_point(size = 3) +
  xlab('Half life, hours') + ylab('Regulated concentration change')+
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),  # Centered, bold title
    axis.text = element_text(color = "black"),  # Black axis text
    axis.title = element_text(face = "bold")  # Bold axis labels
  )+theme_bw(base_size = 15)+
  scale_x_log10(labels = trans_format("log10", math_format(10^.x))) + 
  scale_y_log10(
    labels = trans_format("log10", math_format(10^.x)), 
    breaks = scales::trans_breaks("log10", function(x) 10^x)  # Custom breaks
  ) +
  scale_color_manual(values = c('black','grey50'))+
  scale_alpha_manual(values = c(1,.4))


cor(df_plt$regulated, log2(df_plt$HL),use = 'pairwise.complete.obs')^2/(.60*.95)


##################################################################################################
# Figure 2 e 
##################################################################################################

#Go_Human <- read.delim('/Users/andrewleduc/Desktop/Projects/Miceotopes/Bulk/Gene_sets/GOA_db.txt',sep = ' ')
Go_Human <- read.delim('https://zenodo.org/records/14827610/files/GO_Human.txt?download=1',sep = ' ')

df_plt$vect_norm <- df_plt$vect - median(df_plt$vect,na.rm=T)
df_plt$dilute_expect <- df_plt$V2 - median(df_plt$V2,na.rm=T)


term <- c()
fract_var <- c()
p_val_store <- c()
fc_store <- c()


for(i in unique(Go_Human$GO_term_name)){
  go_hold <- Go_Human %>% filter(GO_term_name == i)
  #go_hold$gene <- toupper(go_hold$gene)
  if(length(intersect(go_hold$Gene,df_plt$gene)) > 3){
    
    df_plt_hold <- df_plt %>% filter(gene %in% go_hold$Gene)
    
    p_val_store <- c(p_val_store,t.test(df_plt_hold$vect_norm,df_plt$vect_norm)$p.value)
    
    
    df_plt_hold <- df_plt_hold %>% filter(is.na(vect_norm)==F)
    df_plt_hold <- df_plt_hold %>% filter(is.na(dilute_expect)==F)

    fract_var <- c(fract_var,median(df_plt_hold$dilute_expect/df_plt_hold$vect_norm,na.rm=T))
    term <- c(term,i)
    fc_store <- c(fc_store,median(df_plt_hold$vect_norm,na.rm=T))
    
  }
}

df_go_bcell <- data.frame(pval =p_val_store,GO = term,FC_avg = fc_store,fract_exp = fract_var )
df_go_bcell$Qval <- p.adjust(df_go_bcell$pval,method = 'BH')

df_go_bcell <- df_go_bcell %>% filter(Qval < .05)
df_go_bcell_ <- df_go_bcell %>% filter(fract_exp > .1)
#write.csv(df_go_bcell_,'~/Desktop/Projects/Dilution_paper/GO_table1.csv')

df_go_bcell$fract_exp_scaled <- (df_go_bcell$fract_exp - min(df_go_bcell$fract_exp, na.rm = TRUE)) /
  (max(df_go_bcell$fract_exp, na.rm = TRUE) - min(df_go_bcell$fract_exp, na.rm = TRUE))

hist(df_go_bcell$fract_exp_scaled)

go_hold <- Go %>% filter(term == 'WP_TCA_CYCLE')
df_plt_hold <- df_plt %>% filter(gene %in% toupper(go_hold$gene))
df_plt_hold$gene
intersect(s,df_plt_hold$gene)

df_plt$col2 <- 'NA'
df_plt$col2[df_plt$gene %in% intersect(df_plt_hold$gene,df_plt$gene)] <- ' AA cat'
#df_plt$col2[df_plt$gene =='USP10'] <- 'IRF4'
ggplot(df_plt,aes(y = 2^(vect - median(vect,na.rm=T)),x = 2^(V2- median(V2,na.rm=T)),color = col2,alpha = col2)) + geom_point(size = 3)+ 
  xlab('Expected dilution') + ylab('Concentration change') +
  geom_abline(intercept = 0, slope = 1)+
  dot_plot + scale_x_log10(labels = trans_format("log10", math_format(10^.x))) + 
  scale_y_log10(labels = trans_format("log10", math_format(10^.x))) +
  scale_color_manual(values = c('black','grey50'))+
  scale_alpha_manual(values = c(1,.1))+
  ggtitle('Branched Chain AA Catabolism')

df_plt$dif_look <- df_plt$vect_norm - df_plt$dilute_expect
df_plt2 <- df_plt %>% filter(abs(dif_look) < .5 )


save_list <- c('HEVNER_CORTEX_PROLIFERATING_CELLS','WP_TRANSLATION_FACTORS',
               'REACTOME_MISMATCH_REPAIR',
               'NABA_ECM_GLYCOPROTEINS','BIOCARTA_NOS1_PATHWAY', 
               'REACTOME_RESPIRATORY_ELECTRON_TRANSPORT',
               'REACTOME_NOTCH_HLH_TRANSCRIPTION_PATHWAY')

df_plt_D <- df_plt %>% dplyr::select(vect_norm,dilute_expect,gene)
View(df_plt)

Go_Terms_boxplot(save_list,df_plt_D)


##################################################################################################
# Figure 2 f
##################################################################################################


df_plt$vect_norm <- df_plt$vect - median(df_plt$vect,na.rm=T)
df_plt$dilute_expect <- df_plt$V2 - median(df_plt$V2,na.rm=T)

df_plt$regulated <- df_plt$vect_norm - df_plt$dilute_expect 

term <- c()
p_val_store <- c()
fc_store <- c()

df_plt <- df_plt %>% filter(is.na(regulated)==F)
for(i in unique(Go$term)){
  go_hold <- Go %>% filter(term == i)
  go_hold$gene <- toupper(go_hold$gene)
  if(length(intersect(go_hold$gene,df_plt$gene)) > 3){
    
    df_plt_hold <- df_plt %>% filter(gene %in% go_hold$gene)
    
    p_val_store <- c(p_val_store,t.test(df_plt_hold$regulated,df_plt$regulated)$p.value)

    df_plt_hold <- df_plt_hold %>% filter(is.na(regulated)==F)

    term <- c(term,i)
    fc_store <- c(fc_store,median(df_plt_hold$regulated - df_plt_hold$vect_norm,na.rm=T))
    
  }
}

df_go_bcell_regulated <- data.frame(pval =p_val_store,GO = term,FC_avg = fc_store )
df_go_bcell_regulated$Qval <- p.adjust(df_go_bcell_regulated$pval,method = 'BH')

df_go_bcell_regulated <- df_go_bcell_regulated %>% filter(Qval < .01)

#write.csv(df_go_bcell_regulated,'~/Desktop/Projects/Dilution_paper/GO_table2.csv')

df_go_bcell_regulated <- df_go_bcell_regulated %>% filter(abs(FC_avg) > .4)
View(df_go_bcell_regulated)

df_plt_plt <- df_plt %>% dplyr::select(vect_norm,regulated,gene)


save_list_reg <- c('REACTOME_GLUCONEOGENESIS','REACTOME_CITRIC_ACID_CYCLE_TCA_CYCLE',
                   'WP_AMINO_ACID_METABOLISM','ISHIDA_E2F_TARGETS')


Go_Terms_boxplot(save_list_reg,df_plt_plt)






#--------------------------------------------------------------------------------------------

# test <- read.delim('/Users/andrewleduc/Library/CloudStorage/GoogleDrive-research@slavovlab.net/My Drive/MS/Users/aleduc/forGeorg/combined/txt/evidence.txt')
# test <- test %>% filter(Acetyl..Protein.N.term. == 0)
# sum(test$TMTPro_Nter_LE==1)/nrow(test)
# length(unique((test$Leading.razor.protein)))
# 
# term <- c()
# p_val_store <- c()
# fc_store <- c()
# 
# count = 1
# for(i in unique(complexs$ComplexName)){
#   go_hold <- complexs %>% filter(ComplexName == i)
#   genes <- unlist(str_split(complexs$subunits.Gene.name.[count],';'))
#   if(length(intersect(genes,df_plt$gene)) > 3){
#     
#     df_plt_hold <- df_plt %>% filter(gene %in% genes)
#     
#     p_val_store <- c(p_val_store,t.test(df_plt_hold$regulated,df_plt_hold$vect_norm)$p.value)
#     
#     df_plt_hold <- df_plt_hold %>% filter(is.na(regulated)==F)
#     
#     term <- c(term,i)
#     fc_store <- c(fc_store,median(df_plt_hold$regulated - df_plt_hold$vect_norm,na.rm=T))
#     
#   }
#   count = 1+count
# }
# 
# df_go_bcell_regulated_complex <- data.frame(pval =p_val_store,GO = term,FC_avg = fc_store )
# df_go_bcell_regulated_complex$Qval <- p.adjust(df_go_bcell_regulated_complex$pval,method = 'BH')
# 
# df_go_bcell_regulated_complex <- df_go_bcell_regulated_complex %>% filter(Qval < .01)
# 
# df_go_bcell_regulated <- df_go_bcell_regulated %>% filter(abs(FC_avg) > .4)
# View(df_go_bcell_regulated)