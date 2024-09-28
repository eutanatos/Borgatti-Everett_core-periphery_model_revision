########################################################################################################################
## REVISING THE BORGATTI-EVERETT CORE-PERIPHERY MODEL
## (4) Simulations
## R script written by José Luis Estévez (University of Helsinki / Vaestoliitto)
## Date: Sep 27th, 2024
########################################################################################################################

# R PACKAGES REQUIRED
library(data.table);library(tidyverse);library(stringr);library(igraph);library(netUtils);library(wCorr);library(lme4)
library(patchwork);library(microbenchmark)

# CLEAN ENVIRONMENT
rm(list=ls())

# LOAD FUNCTIONS
source('03_Functions.R')

# Set theme for plots
theme_set(theme_bw())

########################################################################################################################

# FIRST EXPERIMENT

# Combination of parameters to be check
n_size <- 50 # network size
cvals <- rep(seq(5,45,by=5))
pvals <- 1/9
kvals <- rep(seq(1,3,by=.1))
repetitions <- 100 # number of repetitions per set of parameters

data <- data.table(N = n_size, # number of nodes in the network
                   c = rep(cvals,each=repetitions*length(kvals)), # number of core members
                   p = pvals, # prob of tie, irrespective of block
                   k = rep(kvals,times=repetitions*length(cvals))) # k value

# Synthetic networks
simntw <- list()
set.seed(0708)
for(i in 1:nrow(data)){
  simntw[[i]] <- random_cp(N=data$N[i],c=data$c[i],p=data$p[i],k=data$k[i])
}

# Different implementations applied to synthetic networks
results <- lapply(simntw,cp_ucinet) # standard BE method
results2 <- lapply(simntw,cp_minden,delta=0) # min-density blocks
results3 <- lapply(simntw,cp_pcore,delta=0,p=0.25) # min-density blocks, plus p-core 0.25

# Data extraction
for(i in seq_along(results)){
  # Order of values
  # True negative (predicted and real periphery), 
  # False positive (predicted core, but real periphery)
  # False negative (predicted periphery, but real core)
  # True positive (predicted and real core)
  # Solution and precision (Ucinet style)
  data$ucinet[i] <- paste(as.vector(table(results[[i]]$vec,
                                          as.integer(startsWith(V(simntw[[i]])$name,'C')))),collapse=';')
  # Minimum density blocks
  data$minden[i] <- paste(as.vector(table(results2[[i]]$vec,
                                          as.integer(startsWith(V(simntw[[i]])$name,'C')))),collapse=';')
  # Minimum density blocks and pcore 0.25
  data$pcore[i] <- paste(as.vector(table(results3[[i]]$vec,
                                         as.integer(startsWith(V(simntw[[i]])$name,'C')))),collapse=';')
}

# Long format
data <- data.table(pivot_longer(data,cols=c(ucinet,minden,pcore)))
# Get the values in different columns
vals <- data.table(do.call(rbind,strsplit(data$value,split=';')))
names(vals) <- c('TN','FP','FN','TP') # col names
vals <- apply(vals,2,as.numeric) # values as numeric
data <- cbind(data,vals)

# Save raw data as csv file
write_csv(data,file='simulations1.csv')

# Data wrangling for visualization
mn1 <- data[,mean(TN),by=.(name,k,c)]
mn2 <- data[,mean(FP),by=.(name,k,c)]
mn3 <- data[,mean(FN),by=.(name,k,c)]
mn4 <- data[,mean(TP),by=.(name,k,c)]
# Put altogether
mn1$V2 <- mn2$V1
mn1$V3 <- mn3$V1
mn1$V4 <- mn4$V1
names(mn1) <- c('model','k','c','tn','fp','fn','tp')
# Long format
data <- data.table(pivot_longer(as_tibble(mn1),cols=c(tn,fp,fn,tp)))

# Rename
data[,model := factor(model,levels=c('ucinet','minden','pcore'),
                      labels=c('italic(d)=="NA"','italic(d)>=0','italic(d)>=0 ~~~ italic(p)>=0.25'))]
data[,core_size := factor(c,labels=paste('italic(c)==',cvals,sep=''))]
data[,acc := factor(name,levels=c('fp','tn','tp','fn'),
                    labels=c('False core','True periphery','True core','False periphery'))]

# Visualization
p1.1 <- ggplot(data=data[model %in% c('italic(d)=="NA"','italic(d)>=0')],
               aes(x=k,y=value,fill=acc)) +
  geom_bar(stat='identity') +
  geom_hline(aes(yintercept = c)) +
  geom_text(aes(x=2,y=c+2.5,label='periphery'),size=2.5) +
  geom_text(aes(x=2,y=c-2.5,label='core'),size=2.5) +
  facet_grid(model~core_size,labeller = label_parsed) +
  labs(y='Average confusion matrix',color='',fill='',linetype='') +
  scale_fill_manual(values=c('red','green2','chartreuse','red3')) +
  scale_x_continuous(breaks = seq(1, 3, by = 1)) +
  theme(legend.position="top",
        axis.title.x = element_blank(),   # Remove x-axis label
        axis.text.x = element_blank(),    # Remove x-axis text (tick labels)
        axis.ticks.x = element_blank())

p1.2 <- ggplot(data=data[model %in% c('italic(d)=="NA"','italic(d)>=0 ~~~ italic(p)>=0.25')],
               aes(x=k,y=value,fill=acc)) +
  geom_bar(stat='identity') +
  geom_hline(aes(yintercept = c)) +
  geom_text(aes(x=2,y=c+2.5,label='periphery'),size=2.5) +
  geom_text(aes(x=2,y=c-2.5,label='core'),size=2.5) +
  facet_grid(model~core_size,labeller = label_parsed) +
  labs(y='Average confusion matrix',color='',fill='',linetype='') +
  scale_fill_manual(values=c('red','green2','chartreuse','red3')) +
  scale_x_continuous(breaks = seq(1, 3, by = 1)) +
  theme(legend.position="top",
        axis.title.x = element_blank(),   # Remove x-axis label
        axis.text.x = element_blank(),    # Remove x-axis text (tick labels)
        axis.ticks.x = element_blank())

# Second part of the plot comparing the accuracy of both models
# 2.1) Comparing BE with min-density blocks
data <- data.table(read_csv('simulations1.csv')) # bring data back
# Change ref category
data <- data[name %in% c('ucinet','minden')]
data[,name := factor(name,levels=c('ucinet','minden'))]
# Success proportion
data[,success := (TP+TN)/50]

# Data template
diff.data <- data.table(expand.grid(core_size=unique(data$c),k=unique(data$k)))

for (i in 1:nrow(diff.data)) {
  # Subset of the data with that combination of parameters
  tomodel <- data[c == diff.data$core_size[i] & k == diff.data$k[i]]
  
  # Dynamically create trial_id based on the number of rows in tomodel
  trial_id_count <- nrow(tomodel) / 2
  if (trial_id_count == floor(trial_id_count)) {
    tomodel[, trial_id := rep(1:trial_id_count, each = 2)]
  } else {
    # If trial_id cannot be assigned because of incorrect row counts, skip this iteration
    next
  }
  
  # Try to fit the model
  result <- tryCatch({
    model <- lmer(success ~ -1 + name + (1 | trial_id), data = tomodel)
    
    # Extract estimates and SE
    est <- fixef(model)[1]
    se <- sqrt(diag(vcov(model)))[1]
    est2 <- fixef(model)[2]
    se2 <- sqrt(diag(vcov(model)))[2]
    
    # Return both est and se in a list
    list(est = est, se = se, est2 = est2, se2 = se2)
  }, warning = function(w) {
    # If a warning occurs, return NULL to skip this iteration
    return(NULL)
  }, error = function(e) {
    # If an error occurs, also return NULL
    return(NULL)
  })
  
  # If result is not NULL, store the values in diff.data
  if (!is.null(result)) {
    diff.data$ucinet[i] <- result$est
    diff.data$ucinet_se[i] <- result$se
    diff.data$new[i] <- result$est2
    diff.data$new_se[i] <- result$se2
  } else {
    # Optionally, you can set the values to NA if the model fails
    diff.data$ucinet[i] <- NA
    diff.data$ucinet_se[i] <- NA
    diff.data$new[i] <- NA
    diff.data$new_se[i] <- NA
  }
}

# Change data format
est <- gather(diff.data[,c(1:3,5)],model,est,c(ucinet,new))
se <- gather(diff.data[,c(1:2,4,6)],model,se,c(ucinet_se,new_se))
diff.data <- data.table(est)
diff.data[,se := se$se]
diff.data[,model := factor(model,levels=c('ucinet','new'))]

# Visualization
p2.1 <- ggplot(data=diff.data,aes(x=k,y=est,color=model,fill=model,
                          ymin = est+qnorm(0.025)*se,
                          ymax = est+qnorm(0.975)*se)) +
  geom_ribbon(alpha=.5,linewidth=.1) + 
  geom_line(linewidth=.25) +
  scale_color_manual(values = c('red3','navyblue'), 
                     labels = c(expression(italic(d) == "NA"), expression(italic(d) >= 0))) +
  scale_fill_manual(values = c('orange', 'royalblue'), 
                    labels = c(expression(italic(d) == "NA"), expression(italic(d) >= 0))) +
  facet_grid('Proportion'~core_size,labeller = label_parsed) +
  labs(x=expression(italic(k)),y='Actor assignment accuracy',color='',fill='') +
  scale_x_continuous(breaks = seq(1, 3, by = 1)) +
  theme(strip.text.x = element_blank(),legend.position = 'top')

# Combined both plots
combined_plot <- p1.1 + p2.1 + plot_layout(ncol = 1, heights = c(2, 1))

# Print the combined plot
tiff(filename="Fig7.1.tiff",
     width=27, height=20,units="cm", 
     compression="lzw",bg="white",res=1000
)
combined_plot
dev.off()

# 2.2) Comparing BE with combiantion min-density and p-core 0.25
data <- data.table(read_csv('simulations1.csv')) # bring data back
# Change ref category
data <- data[name %in% c('ucinet','pcore')]
data[,name := factor(name,levels=c('ucinet','pcore'))]
# Success proportion
data[,success := (TP+TN)/50]

# Data template
diff.data <- data.table(expand.grid(core_size=unique(data$c),k=unique(data$k)))

for (i in 1:nrow(diff.data)) {
  # Subset of the data with that combination of parameters
  tomodel <- data[c == diff.data$core_size[i] & k == diff.data$k[i]]
  
  # Dynamically create trial_id based on the number of rows in tomodel
  trial_id_count <- nrow(tomodel) / 2
  if (trial_id_count == floor(trial_id_count)) {
    tomodel[, trial_id := rep(1:trial_id_count, each = 2)]
  } else {
    # If trial_id cannot be assigned because of incorrect row counts, skip this iteration
    next
  }
  
  # Try to fit the model
  result <- tryCatch({
    model <- lmer(success ~ -1 + name + (1 | trial_id), data = tomodel)
    
    # Extract estimates and SE
    est <- fixef(model)[1]
    se <- sqrt(diag(vcov(model)))[1]
    est2 <- fixef(model)[2]
    se2 <- sqrt(diag(vcov(model)))[2]
    
    # Return both est and se in a list
    list(est = est, se = se, est2 = est2, se2 = se2)
  }, warning = function(w) {
    # If a warning occurs, return NULL to skip this iteration
    return(NULL)
  }, error = function(e) {
    # If an error occurs, also return NULL
    return(NULL)
  })
  
  # If result is not NULL, store the values in diff.data
  if (!is.null(result)) {
    diff.data$ucinet[i] <- result$est
    diff.data$ucinet_se[i] <- result$se
    diff.data$new[i] <- result$est2
    diff.data$new_se[i] <- result$se2
  } else {
    # Optionally, you can set the values to NA if the model fails
    diff.data$ucinet[i] <- NA
    diff.data$ucinet_se[i] <- NA
    diff.data$new[i] <- NA
    diff.data$new_se[i] <- NA
  }
}

# Change data format
est <- gather(diff.data[,c(1:3,5)],model,est,c(ucinet,new))
se <- gather(diff.data[,c(1:2,4,6)],model,se,c(ucinet_se,new_se))
diff.data <- data.table(est)
diff.data[,se := se$se]
diff.data[,model := factor(model,levels=c('ucinet','new'))]

# Visualization
p2.2 <- ggplot(data=diff.data,aes(x=k,y=est,color=model,fill=model,
                          ymin = est+qnorm(0.025)*se,
                          ymax = est+qnorm(0.975)*se)) +
  geom_ribbon(alpha=.5,linewidth=.1) + 
  geom_line(linewidth=.25) +
  scale_color_manual(values = c('red3','navyblue'), 
                     labels = c(expression(italic(d) == "NA"), expression(italic(d) >= 0 ~ "&" ~ italic(p) >= 0.25))) +
  scale_fill_manual(values = c('orange', 'royalblue'), 
                    labels = c(expression(italic(d) == "NA"), expression(italic(d) >= 0 ~ "&" ~ italic(p) >= 0.25))) +
  facet_grid('Proportion'~core_size,labeller = label_parsed) +
  labs(x=expression(italic(k)),y='Actor assignment accuracy',color='',fill='') +
  scale_x_continuous(breaks = seq(1, 3, by = 1)) +
  theme(strip.text.x = element_blank(),legend.position = 'top')

# Combined both plots
combined_plot2 <- p1.2 + p2.2 + plot_layout(ncol = 1, heights = c(2, 1))

# Print the combined plot
tiff(filename="Fig7.2.tiff",
     width=27, height=20,units="cm", 
     compression="lzw",bg="white",res=1000
)
combined_plot2
dev.off()

###########################

# Comparison of computational costs

# Use a random sample (100) of the 9,450 synthetic networks
set.seed(0708)
sample <- sample(1:length(simntw),100,replace=FALSE)

result1 <- microbenchmark(
  lapply(simntw[sample],cp_ucinet),
  times = 100  # Number of times to repeat the measurement
)
print(result1)

result2 <- microbenchmark(
  lapply(simntw[sample],cp_minden,delta=0),
  times = 100  # Number of times to repeat the measurement
)
print(result2)

result3 <- microbenchmark(
  lapply(simntw[sample],cp_pcore,delta=0,p=0.25),
  times = 100  # Number of times to repeat the measurement
)
print(result3)

########################################################################################################################

# SECOND EXPERIMENT

# Illustrative examples
set.seed(0708)
proxies <- list()
proxies[[1]] <- random_cp2(N=50,c=20,p=1/9,k=3,connected='no')
proxies[[2]] <- random_cp2(N=50,c=20,p=1/9,k=3,connected='core')
proxies[[3]] <- random_cp2(N=50,c=20,p=1/9,k=3,connected='periphery')
proxies[[4]] <- random_cp2(N=50,c=20,p=1/9,k=2,connected='no')
proxies[[5]] <- random_cp2(N=50,c=20,p=1/9,k=2,connected='core')
proxies[[6]] <- random_cp2(N=50,c=20,p=1/9,k=2,connected='periphery')

titles <- list(
  expression(paste("Unconnected, ", italic(c) == 20,', ',italic(k) == 3)),
  expression(paste("Core-connected, ", italic(c) == 20,', ',italic(k) == 3)),
  expression(paste("Periphery-connected, ", italic(c) == 20,', ',italic(k) == 3)),
  expression(paste("Unconnected, ", italic(c) == 20,', ',italic(k) == 2)),
  expression(paste("Core-connected, ", italic(c) == 20,', ',italic(k) == 2)),
  expression(paste("Periphery-connected, ", italic(c) == 20,', ',italic(k) == 2))
)

tiff(filename="Fig8.tiff",
     width=20, height=15,units="cm", 
     compression="lzw",bg="white",res=1000
)
par(mfrow=c(2,3))
for(i in seq_along(proxies)){
  plot.igraph(proxies[[i]],
              vertex.label = '',
              vertex.color=ifelse(startsWith(V(proxies[[1]])$name,'C'),'grey10','grey95'),
              edge.color='grey70',edge.width=1.5,
              layout = layout_with_kk(proxies[[i]]),
              main = titles[i])
}
dev.off()

##############################

# Combination of parameters to be check
n_size <- 50 # network size
cvals <- c(10,20,30,40)
pvals <- 1/9
kvals <- rep(seq(1,3,by=.1))
conn <- c('no','core','periphery')
repetitions <- 100 # number of repetitions per set of parameters

data <- data.table(N = n_size, # number of nodes in the network
                   c = rep(cvals,each=repetitions*length(kvals)*length(conn)), # number of core members
                   p = pvals, # prob of tie, irrespective of block
                   k = rep(kvals,times=repetitions*length(cvals)*length(conn)), # k value
                   conn = rep(conn,times=repetitions*length(cvals),each=length(kvals))) # type of connectivity

# Synthetic networks
simntw <- list()
set.seed(0708)
for(i in 1:nrow(data)){
  simntw[[i]] <- random_cp2(N=data$N[i],c=data$c[i],p=data$p[i],k=data$k[i],connected=data$conn[i])
}

# Different implementations applied to synthetic networks
results <- lapply(simntw,cp_ucinet) # standard BE method
results2 <- lapply(simntw,cp_minden,delta=0) # min-density blocks
results3 <- lapply(simntw,cp_pcore,delta=0,p=0.25) # min-density blocks, plus p-core 0.25

# Data extraction
for(i in seq_along(results)){
  # Solution and precision (Ucinet style)
  data$ucinet[i] <- paste(as.vector(table(results[[i]]$vec,
                                          as.integer(startsWith(V(simntw[[i]])$name,'C')))),collapse=';')
  # Minimum density blocks
  data$minden[i] <- paste(as.vector(table(results2[[i]]$vec,
                                          as.integer(startsWith(V(simntw[[i]])$name,'C')))),collapse=';')
  # Minimum density blocks and pcore 0.25
  data$pcore[i] <- paste(as.vector(table(results3[[i]]$vec,
                                         as.integer(startsWith(V(simntw[[i]])$name,'C')))),collapse=';')
}

# Long format
data <- data.table(pivot_longer(data,cols=c(ucinet,minden,pcore)))
# Get the values in different columns
vals <- data.table(do.call(rbind,strsplit(data$value,split=';')))
names(vals) <- c('TN','FP','FN','TP') # col names
vals <- apply(vals,2,as.numeric) # values as numeric
data <- cbind(data,vals)

# Save raw data as csv file
write_csv(data,file='simulations2.csv')

# Data wrangling for visualization
mn1 <- data[,mean(TN),by=.(name,k,c,conn)]
mn2 <- data[,mean(FP),by=.(name,k,c,conn)]
mn3 <- data[,mean(FN),by=.(name,k,c,conn)]
mn4 <- data[,mean(TP),by=.(name,k,c,conn)]
# Put altogether
mn1$V2 <- mn2$V1
mn1$V3 <- mn3$V1
mn1$V4 <- mn4$V1
names(mn1) <- c('model','k','c','conn','tn','fp','fn','tp')
# Long format
data <- data.table(pivot_longer(as_tibble(mn1),cols=c(tn,fp,fn,tp)))

# Rename
data[,model := factor(model,levels=c('ucinet','minden','pcore'),
                      labels=c('italic(d)=="NA"','italic(d)>=0','italic(d)>=0 ~~~ italic(p)>=0.25'))]
data[,core_size := factor(c,labels=paste('italic(c)==',cvals,sep=''))]
data[,acc := factor(name,levels=c('fp','tn','tp','fn'),
                    labels=c('False core','True periphery','True core','False periphery'))]
data[,type := factor(conn,levels=c('no','core','periphery'),labels=c('Unconnected','Core','Periphery'))]

# Visualization
p3.1 <- ggplot(data=data[model %in% c('italic(d)=="NA"','italic(d)>=0')],
               aes(x=k,y=value,fill=acc)) +
  geom_bar(stat='identity') +
  geom_hline(aes(yintercept = c)) +
  geom_text(aes(x=2,y=c+2.5,label='periphery'),size=2.5) +
  geom_text(aes(x=2,y=c-2.5,label='core'),size=2.5) +
  facet_grid(model~core_size+type,labeller = label_parsed) +
  labs(x=expression(italic(k)),y='Average confusion matrix',color='',fill='',linetype='') +
  scale_fill_manual(values=c('red','green2','chartreuse','red3')) +
  scale_x_continuous(breaks = seq(1, 3, by = 1)) +
  theme(legend.position="top",
        axis.title.x = element_blank(),   # Remove x-axis label
        axis.text.x = element_blank(),    # Remove x-axis text (tick labels)
        axis.ticks.x = element_blank())

p3.2 <- ggplot(data=data[model %in% c('italic(d)=="NA"','italic(d)>=0 ~~~ italic(p)>=0.25')],
               aes(x=k,y=value,fill=acc)) +
  geom_bar(stat='identity') +
  geom_hline(aes(yintercept = c)) +
  geom_text(aes(x=2,y=c+2.5,label='periphery'),size=2.5) +
  geom_text(aes(x=2,y=c-2.5,label='core'),size=2.5) +
  facet_grid(model~core_size+type,labeller = label_parsed) +
  labs(x=expression(italic(k)),y='Average confusion matrix',color='',fill='',linetype='') +
  scale_fill_manual(values=c('red','green2','chartreuse','red3')) +
  scale_x_continuous(breaks = seq(1, 3, by = 1)) +
  theme(legend.position="top",
        axis.title.x = element_blank(),   # Remove x-axis label
        axis.text.x = element_blank(),    # Remove x-axis text (tick labels)
        axis.ticks.x = element_blank())

# Second part of the plot comparing the accuracy of both models
# 2.1) Comparing BE with min-density blocks
data <- data.table(read_csv('simulations2.csv')) # bring data back
# Change ref category
data <- data[name %in% c('ucinet','minden')]
data[,name := factor(name,levels=c('ucinet','minden'))]
# Success proportion
data[,success := (TP+TN)/50]

# Data template
diff.data <- data.table(expand.grid(core_size=unique(data$c),k=unique(data$k),conn=unique(data$conn)))

for (i in 1:nrow(diff.data)) {
  # Subset of the data with that combination of parameters
  tomodel <- data[c == diff.data$core_size[i] & k == diff.data$k[i] & conn == diff.data$conn[i]]
  
  # Dynamically create trial_id based on the number of rows in tomodel
  trial_id_count <- nrow(tomodel) / 2
  if (trial_id_count == floor(trial_id_count)) {
    tomodel[, trial_id := rep(1:trial_id_count, each = 2)]
  } else {
    # If trial_id cannot be assigned because of incorrect row counts, skip this iteration
    next
  }
  
  # Try to fit the model
  result <- tryCatch({
    model <- lmer(success ~ -1 + name + (1 | trial_id), data = tomodel)
    
    # Extract estimates and SE
    est <- fixef(model)[1]
    se <- sqrt(diag(vcov(model)))[1]
    est2 <- fixef(model)[2]
    se2 <- sqrt(diag(vcov(model)))[2]
    
    # Return both est and se in a list
    list(est = est, se = se, est2 = est2, se2 = se2)
  }, warning = function(w) {
    # If a warning occurs, return NULL to skip this iteration
    return(NULL)
  }, error = function(e) {
    # If an error occurs, also return NULL
    return(NULL)
  })
  
  # If result is not NULL, store the values in diff.data
  if (!is.null(result)) {
    diff.data$ucinet[i] <- result$est
    diff.data$ucinet_se[i] <- result$se
    diff.data$new[i] <- result$est2
    diff.data$new_se[i] <- result$se2
  } else {
    # Optionally, you can set the values to NA if the model fails
    diff.data$ucinet[i] <- NA
    diff.data$ucinet_se[i] <- NA
    diff.data$new[i] <- NA
    diff.data$new_se[i] <- NA
  }
}

# Change data format
est <- gather(diff.data[,c(1:4,6)],model,est,c(ucinet,new))
se <- gather(diff.data[,c(1:3,5,7)],model,se,c(ucinet_se,new_se))
diff.data <- data.table(est)
diff.data[,se := se$se]
diff.data[,model := factor(model,levels=c('ucinet','new'))]

# Visualization
p4.1 <- ggplot(data=diff.data,aes(x=k,y=est,color=model,fill=model,
                          ymin = est+qnorm(0.025)*se,
                          ymax = est+qnorm(0.975)*se)) +
  geom_ribbon(alpha=.5,linewidth=.1) + 
  geom_line(linewidth=.25) +
  scale_color_manual(values = c('red3','navyblue'), 
                     labels = c(expression(italic(d) == "NA"), expression(italic(d) >= 0))) +
  scale_fill_manual(values = c('orange', 'royalblue'), 
                    labels = c(expression(italic(d) == "NA"), expression(italic(d) >= 0))) +
  facet_grid('Proportion'~core_size+conn,labeller = label_parsed) +
  labs(x=expression(italic(k)),y='Actor assignment accuracy',color='',fill='') +
  scale_x_continuous(breaks = seq(1, 3, by = 1)) +
  theme(strip.text.x = element_blank(),legend.position = 'top')

# Combined both plots
combined_plot3 <- p3.1 + p4.1 + plot_layout(ncol = 1, heights = c(2, 1))

# Print the combined plot
tiff(filename="Fig9.1.tiff",
     width=27, height=20,units="cm", 
     compression="lzw",bg="white",res=1000
)
combined_plot3
dev.off()

# 2.2) Comparing BE with min-density blocks and pcore 0.25
data <- data.table(read_csv('simulations2.csv')) # bring data back
# Change ref category
data <- data[name %in% c('ucinet','pcore')]
data[,name := factor(name,levels=c('ucinet','pcore'))]
# Success proportion
data[,success := (TP+TN)/50]

# Data template
diff.data <- data.table(expand.grid(core_size=unique(data$c),k=unique(data$k),conn=unique(data$conn)))

for (i in 1:nrow(diff.data)) {
  # Subset of the data with that combination of parameters
  tomodel <- data[c == diff.data$core_size[i] & k == diff.data$k[i] & conn == diff.data$conn[i]]
  
  # Dynamically create trial_id based on the number of rows in tomodel
  trial_id_count <- nrow(tomodel) / 2
  if (trial_id_count == floor(trial_id_count)) {
    tomodel[, trial_id := rep(1:trial_id_count, each = 2)]
  } else {
    # If trial_id cannot be assigned because of incorrect row counts, skip this iteration
    next
  }
  
  # Try to fit the model
  result <- tryCatch({
    model <- lmer(success ~ -1 + name + (1 | trial_id), data = tomodel)
    
    # Extract estimates and SE
    est <- fixef(model)[1]
    se <- sqrt(diag(vcov(model)))[1]
    est2 <- fixef(model)[2]
    se2 <- sqrt(diag(vcov(model)))[2]
    
    # Return both est and se in a list
    list(est = est, se = se, est2 = est2, se2 = se2)
  }, warning = function(w) {
    # If a warning occurs, return NULL to skip this iteration
    return(NULL)
  }, error = function(e) {
    # If an error occurs, also return NULL
    return(NULL)
  })
  
  # If result is not NULL, store the values in diff.data
  if (!is.null(result)) {
    diff.data$ucinet[i] <- result$est
    diff.data$ucinet_se[i] <- result$se
    diff.data$new[i] <- result$est2
    diff.data$new_se[i] <- result$se2
  } else {
    # Optionally, you can set the values to NA if the model fails
    diff.data$ucinet[i] <- NA
    diff.data$ucinet_se[i] <- NA
    diff.data$new[i] <- NA
    diff.data$new_se[i] <- NA
  }
}

# Change data format
est <- gather(diff.data[,c(1:4,6)],model,est,c(ucinet,new))
se <- gather(diff.data[,c(1:3,5,7)],model,se,c(ucinet_se,new_se))
diff.data <- data.table(est)
diff.data[,se := se$se]
diff.data[,model := factor(model,levels=c('ucinet','new'))]

# Visualization
p4.2 <- ggplot(data=diff.data,aes(x=k,y=est,color=model,fill=model,
                          ymin = est+qnorm(0.025)*se,
                          ymax = est+qnorm(0.975)*se)) +
  geom_ribbon(alpha=.5,linewidth=.1) + 
  geom_line(linewidth=.25) +
  scale_color_manual(values = c('red3','navyblue'), 
                     labels = c(expression(italic(d) == "NA"), expression(italic(d) >= 0 ~ "&" ~ italic(p) >= 0.25))) +
  scale_fill_manual(values = c('orange', 'royalblue'), 
                    labels = c(expression(italic(d) == "NA"), expression(italic(d) >= 0 ~ "&" ~ italic(p) >= 0.25))) +
  facet_grid('Proportion'~core_size+conn,labeller = label_parsed) +
  labs(x=expression(italic(k)),y='Actor assignment accuracy',color='',fill='') +
  scale_x_continuous(breaks = seq(1, 3, by = 1)) +
  theme(strip.text.x = element_blank(),legend.position = 'top')

# Combined both plots
combined_plot4 <- p3.2 + p4.2 + plot_layout(ncol = 1, heights = c(2, 1))

# Print the combined plot
tiff(filename="Fig9.2.tiff",
     width=27, height=20,units="cm", 
     compression="lzw",bg="white",res=1000
)
combined_plot4
dev.off()

########################################################################################################################