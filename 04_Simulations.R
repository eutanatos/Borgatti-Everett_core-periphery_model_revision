################################################################################
## REVISING THE BORGATTI-EVERETT CORE-PERIPHERY MODEL
## (4) Simulations
## R script written by José Luis Estévez (University of Helsinki)
## Date: Nov 21st, 2024
################################################################################

# R PACKAGES REQUIRED
library(data.table);library(tidyverse);library(stringr);library(patchwork);library(scales)
library(igraph);library(netUtils);library(wCorr)
library(lme4);library(sjPlot)
library(microbenchmark)

# CLEAN ENVIRONMENT
rm(list=ls())

# LOAD FUNCTIONS
source('03_Functions.R')

# Set theme for plots
theme_set(theme_bw())

################################################################################

# FIRST EXPERIMENT ----

# Combination of parameters to be check
n_size <- 50 # network size
cvals <- rep(seq(5,45,by=5))
pvals <- 1/9
kvals <- rep(seq(1,3,by=.1))
repetitions <- 50 # number of repetitions per set of parameters

data <- data.table(N = n_size, # number of nodes in the network
                   c = rep(cvals,each=repetitions*length(kvals)), # number of core members
                   p = pvals, # prob of tie, irrespective of block
                   k = rep(kvals,times=repetitions*length(cvals))) # k value

# Synthetic networks
simntw <- list()
set.seed(0708)
for(i in 1:nrow(data)){
  simntw[[i]] <- random.cp(N=data$N[i],c=data$c[i],p=data$p[i],k=data$k[i])
}

# Different implementations applied to synthetic networks
results <- lapply(simntw,cp.ucinet) # standard BE method
results2 <- lapply(simntw,cp.minden,delta=0) # min-density blocks
results3 <- lapply(simntw,cp.pcore,delta=NA,p=0.5) # p-core 0.5

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
  # pcore 0.5
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
                      labels=c('italic(d)=="NA"','italic(d)==0','italic(d)=="NA" ~~~ italic(p)==0.5'))]
data[,core_size := factor(c,labels=paste('italic(c)==',cvals,sep=''))]
data[,acc := factor(name,levels=c('fp','tn','tp','fn'),
                    labels=c('False core','True periphery','True core','False periphery'))]

# Visualization
p1.1 <- ggplot(data=data[model %in% c('italic(d)=="NA"','italic(d)==0')],
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

p1.2 <- ggplot(data=data[model %in% c('italic(d)=="NA"','italic(d)=="NA" ~~~ italic(p)==0.5')],
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
# 2.1) Comparing BE with min-density blocks ----
data <- data.table(read_csv('simulations1.csv')) # bring data back
# Change ref category
data <- data[name %in% c('ucinet','minden')]
data[,name := factor(name,levels=c('ucinet','minden'))]
# Add indicator of network
data[,ntw := rep(1:(nrow(data)/2),each=2)]
data[,ntw := as.factor(ntw)]
# Let's model c as a factor (a set of dummies)
data[,c := as.factor(c)]
# Successes and failures
data[,success := TP+TN]
data[,fail := FP+FN]

# Data modelling
model1 <- glm(cbind(success, fail) ~ name*k*c,
                family = binomial(link = 'logit'), data = data)
# Extract predicted values
plot_model(model1,terms=c('name'),type='pred')
(p2.1 <- plot_model(model1,terms=c('k [all]','name','c'),type='pred'))
# Let's just customize the plots
data <- data.table(p2.1$data)
p2.1 <- ggplot(data=data,aes(x=x,y=predicted,ymin = conf.low, ymax = conf.high,
                     color=group_col,fill=group_col)) + 
  geom_ribbon(alpha=.5,linewidth=.1) + geom_line(linewidth=.25) +
  scale_color_manual(values = c('red3','navyblue'), 
                     labels = c(expression(italic(d) == "NA"), expression(italic(d) == 0))) +
  scale_fill_manual(values = c('orange', 'royalblue'), 
                    labels = c(expression(italic(d) == "NA"), expression(italic(d) == 0 ))) +
  facet_grid('Predicted probabilities'~facet) +
  labs(x=expression(italic(k)),y='Actor assignment accuracy',color='',fill='') +
  scale_y_continuous(labels = percent_format(scale = 100)) +
  scale_x_continuous(breaks = seq(1, 3, by = 1)) +
  theme(strip.text.x = element_blank(),legend.position = 'top')
 
# Combined both plots
combined_plot <- p1.1 + p2.1 + plot_layout(ncol = 1, heights = c(2, 1))

# Print the combined plot
tiff(filename="Fig5.tiff",
     width=27, height=20,units="cm", 
     compression="lzw",bg="white",res=1000
)
combined_plot
dev.off()

# 2.2) Comparing BE with combination min-density and p-core 0.25 ----
data <- data.table(read_csv('simulations1.csv')) # bring data back
# Change ref category
data <- data[name %in% c('ucinet','pcore')]
data[,name := factor(name,levels=c('ucinet','pcore'))]
# Add indicator of network
data[,ntw := rep(1:(nrow(data)/2),each=2)]
data[,ntw := as.factor(ntw)]
# Let's model c as a factor (a set of dummies)
data[,c := as.factor(c)]
# Successes and failures
data[,success := TP+TN]
data[,fail := FP+FN]

# Data modelling
model2 <- glm(cbind(success, fail) ~ name*k*c,
              family = binomial(link = 'logit'), data = data)
# Extract predicted values
plot_model(model2,terms=c('name'),type='pred')
(p2.2 <- plot_model(model2,terms=c('k [all]','name','c'),type='pred'))
# Let's just customize the plots
data <- data.table(p2.2$data)
p2.2 <- ggplot(data=data,aes(x=x,y=predicted,ymin = conf.low, ymax = conf.high,
                             color=group_col,fill=group_col)) + 
  geom_ribbon(alpha=.5,linewidth=.1) + geom_line(linewidth=.25) +
  scale_color_manual(values = c('red3','navyblue'), 
                     labels = c(expression(italic(d) == "NA"), expression(italic(d) == "NA" ~ "&" ~ italic(p) == 0.5))) +
  scale_fill_manual(values = c('orange', 'royalblue'), 
                    labels = c(expression(italic(d) == "NA"), expression(italic(d) == "NA" ~ "&" ~ italic(p) == 0.5))) +
  facet_grid('Predicted probabilities'~facet) +
  labs(x=expression(italic(k)),y='Actor assignment accuracy',color='',fill='') +
  scale_y_continuous(labels = percent_format(scale = 100)) +
  scale_x_continuous(breaks = seq(1, 3, by = 1)) +
  theme(strip.text.x = element_blank(),legend.position = 'top')

# Combined both plots
combined_plot2 <- p1.2 + p2.2 + plot_layout(ncol = 1, heights = c(2, 1))

# Print the combined plot
tiff(filename="Fig8.tiff",
     width=27, height=20,units="cm", 
     compression="lzw",bg="white",res=1000
)
combined_plot2
dev.off()

###########################

# Comparison of computational costs ---

# Use a random sample (100) of the 9,450 synthetic networks
set.seed(0708)
sample <- sample(1:length(simntw),100,replace=FALSE)

result1 <- microbenchmark(
  lapply(simntw[sample],cp.ucinet),
  times = 100  # Number of times to repeat the measurement
)
print(result1)

result2 <- microbenchmark(
  lapply(simntw[sample],cp.minden,delta=0),
  times = 100  # Number of times to repeat the measurement
)
print(result2)

result3 <- microbenchmark(
  lapply(simntw[sample],cp.pcore,delta=NA,p=0.5),
  times = 100  # Number of times to repeat the measurement
)
print(result3)

################################################################################

# SECOND EXPERIMENT ----

# Illustrative examples
set.seed(0708)
proxies <- list()
proxies[[1]] <- random.cp2(N=50,c=20,p=1/9,k=3,connected='no')
proxies[[2]] <- random.cp2(N=50,c=20,p=1/9,k=3,connected='core')
proxies[[3]] <- random.cp2(N=50,c=20,p=1/9,k=3,connected='periphery')
proxies[[4]] <- random.cp2(N=50,c=20,p=1/9,k=2,connected='no')
proxies[[5]] <- random.cp2(N=50,c=20,p=1/9,k=2,connected='core')
proxies[[6]] <- random.cp2(N=50,c=20,p=1/9,k=2,connected='periphery')

titles <- list(
  expression(paste("Unconnected, ", italic(c) == 20,', ',italic(k) == 3)),
  expression(paste("Core-connected, ", italic(c) == 20,', ',italic(k) == 3)),
  expression(paste("Periphery-connected, ", italic(c) == 20,', ',italic(k) == 3)),
  expression(paste("Unconnected, ", italic(c) == 20,', ',italic(k) == 2)),
  expression(paste("Core-connected, ", italic(c) == 20,', ',italic(k) == 2)),
  expression(paste("Periphery-connected, ", italic(c) == 20,', ',italic(k) == 2))
)

tiff(filename="FigA2.tiff",
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
repetitions <- 50 # number of repetitions per set of parameters

data <- data.table(N = n_size, # number of nodes in the network
                   c = rep(cvals,each=repetitions*length(kvals)*length(conn)), # number of core members
                   p = pvals, # prob of tie, irrespective of block
                   k = rep(kvals,times=repetitions*length(cvals)*length(conn)), # k value
                   conn = rep(conn,times=repetitions*length(cvals),each=length(kvals))) # type of connectivity

# Synthetic networks
simntw <- list()
set.seed(0708)
for(i in 1:nrow(data)){
  simntw[[i]] <- random.cp2(N=data$N[i],c=data$c[i],p=data$p[i],k=data$k[i],connected=data$conn[i])
}

# Different implementations applied to synthetic networks
results <- lapply(simntw,cp.ucinet) # standard BE method
results2 <- lapply(simntw,cp.minden,delta=0) # min-density blocks
results3 <- lapply(simntw,cp.pcore,delta=NA,p=0.5) # p-core 0.5

# Data extraction
for(i in seq_along(results)){
  # Solution and precision (Ucinet style)
  data$ucinet[i] <- paste(as.vector(table(results[[i]]$vec,
                                          as.integer(startsWith(V(simntw[[i]])$name,'C')))),collapse=';')
  # Minimum density blocks
  data$minden[i] <- paste(as.vector(table(results2[[i]]$vec,
                                          as.integer(startsWith(V(simntw[[i]])$name,'C')))),collapse=';')
  # p-core 0.5
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
                      labels=c('italic(d)=="NA"','italic(d)==0','italic(d)=="NA" ~~~ italic(p)==0.5'))]
data[,core_size := factor(c,labels=paste('italic(c)==',cvals,sep=''))]
data[,acc := factor(name,levels=c('fp','tn','tp','fn'),
                    labels=c('False core','True periphery','True core','False periphery'))]
data[,type := factor(conn,levels=c('no','core','periphery'),labels=c('Unconnected','Core','Periphery'))]

# Visualization
p3.1 <- ggplot(data=data[model %in% c('italic(d)=="NA"','italic(d)==0')],
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

p3.2 <- ggplot(data=data[model %in% c('italic(d)=="NA"','italic(d)=="NA" ~~~ italic(p)==0.5')],
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
# 2.1) Comparing BE with min-density blocks ----
data <- data.table(read_csv('simulations2.csv')) # bring data back
# Change ref category
data <- data[name %in% c('ucinet','minden')]
data[,name := factor(name,levels=c('ucinet','minden'))]
# Add indicator of network
data[,ntw := rep(1:(nrow(data)/2),each=2)]
data[,ntw := as.factor(ntw)]
# Let's model c as a factor (a set of dummies)
data[,c := as.factor(c)]
# Also redefine the levels of conn
data[,conn := factor(conn,levels=c('no','core','periphery'))]
# Successes and failures
data[,success := TP+TN]
data[,fail := FP+FN]

# Data modelling
model3 <- glm(cbind(success, fail) ~ name*k*c*conn,
              family = binomial(link = 'logit'), data = data)
# Extract predicted values
plot_model(model3,terms=c('name'),type='pred')
(p4.1 <- plot_model(model3,terms=c('k [all]','name','conn','c'),type='pred'))
# Let's just customize the plots
data <- rbind(data.table(p4.1[[1]]$data),data.table(p4.1[[2]]$data),
              data.table(p4.1[[3]]$data),data.table(p4.1[[4]]$data))
p4.1 <- ggplot(data=data,aes(x=x,y=predicted,ymin = conf.low, ymax = conf.high,
                             color=group_col,fill=group_col)) + 
  geom_ribbon(alpha=.5,linewidth=.1) + geom_line(linewidth=.25) +
  scale_color_manual(values = c('red3','navyblue'), 
                     labels = c(expression(italic(d) == "NA"), expression(italic(d) == 0))) +
  scale_fill_manual(values = c('orange', 'royalblue'), 
                    labels = c(expression(italic(d) == "NA"), expression(italic(d) == 0 ))) +
  facet_grid('Predicted probabilities'~panel+facet) +
  labs(x=expression(italic(k)),y='Actor assignment accuracy',color='',fill='') +
  scale_y_continuous(labels = percent_format(scale = 100)) +
  scale_x_continuous(breaks = seq(1, 3, by = 1)) +
  theme(strip.text.x = element_blank(),legend.position = 'top')

# Combined both plots
combined_plot3 <- p3.1 + p4.1 + plot_layout(ncol = 1, heights = c(2, 1))

# Print the combined plot
tiff(filename="FigA3.tiff",
     width=27, height=20,units="cm", 
     compression="lzw",bg="white",res=1000
)
combined_plot3
dev.off()

# 2.2) Comparing BE with min-density blocks and pcore 0.5 ----
data <- data.table(read_csv('simulations2.csv')) # bring data back
# Change ref category
data <- data[name %in% c('ucinet','pcore')]
data[,name := factor(name,levels=c('ucinet','pcore'))]
# Add indicator of network
data[,ntw := rep(1:(nrow(data)/2),each=2)]
data[,ntw := as.factor(ntw)]
# Let's model c as a factor (a set of dummies)
data[,c := as.factor(c)]
# Also redefine the levels of conn
data[,conn := factor(conn,levels=c('no','core','periphery'))]
# Successes and failures
data[,success := TP+TN]
data[,fail := FP+FN]

# Data modelling
model4 <- glm(cbind(success, fail) ~ name*k*c*conn,
              family = binomial(link = 'logit'), data = data)
# Extract predicted values
plot_model(model4,terms=c('name'),type='pred')
(p4.2 <- plot_model(model4,terms=c('k [all]','name','conn','c'),type='pred'))
# Let's just customize the plots
data <- rbind(data.table(p4.2[[1]]$data),data.table(p4.2[[2]]$data),
              data.table(p4.2[[3]]$data),data.table(p4.2[[4]]$data))
p4.2 <- ggplot(data=data,aes(x=x,y=predicted,ymin = conf.low, ymax = conf.high,
                             color=group_col,fill=group_col)) + 
  geom_ribbon(alpha=.5,linewidth=.1) + geom_line(linewidth=.25) +
  scale_color_manual(values = c('red3','navyblue'), 
                     labels = c(expression(italic(d) == "NA"), expression(italic(d) == "NA" ~ "&" ~ italic(p) == 0.5))) +
  scale_fill_manual(values = c('orange', 'royalblue'), 
                    labels = c(expression(italic(d) == "NA"), expression(italic(d) == "NA" ~ "&" ~ italic(p) == 0.5))) +
  facet_grid('Predicted probabilities'~panel+facet) +
  labs(x=expression(italic(k)),y='Actor assignment accuracy',color='',fill='') +
  scale_y_continuous(labels = percent_format(scale = 100)) +
  scale_x_continuous(breaks = seq(1, 3, by = 1)) +
  theme(strip.text.x = element_blank(),legend.position = 'top')

# Combined both plots
combined_plot4 <- p3.2 + p4.2 + plot_layout(ncol = 1, heights = c(2, 1))

# Print the combined plot
tiff(filename="FigA4.tiff",
     width=27, height=20,units="cm", 
     compression="lzw",bg="white",res=1000
)
combined_plot4
dev.off()

################################################################################