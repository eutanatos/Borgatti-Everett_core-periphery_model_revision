################################################################################
## REVISING THE BORGATTI-EVERETT CORE-PERIPHERY MODEL
## (4) Simulations
## R script written by José Luis Estévez (University of Helsinki / Vaestoliitto)
## Date: Oct 11th, 2024
################################################################################

# R PACKAGES REQUIRED
library(data.table);library(tidyverse);library(stringr);library(patchwork);library(scales)
library(igraph);library(netUtils);library(wCorr)
library(glmmTMB);library(ggeffects);library(performance)
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
results3 <- lapply(simntw,cp.pcore,delta=0,p=0.25) # min-density blocks, plus p-core 0.25

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

# Data wrangling for visualization ----
means <- data[, lapply(.SD, mean), .SDcols = c("TN", "FP", "FN", "TP"), by = .(name, k, c)]
setnames(means, c("TN", "FP", "FN", "TP"), paste0("mean_", c("tn", "fp", "fn", "tp")))

# Reshape data to long format
data_long <- as.data.table(pivot_longer(as_tibble(means), cols = starts_with("mean_"), names_to = "acc", values_to = "value"))

# Rename factors
data_long[, `:=`(
  model = factor(name, levels = c("ucinet", "minden", "pcore"),
                 labels = c("italic(d) == 'NA'", 
                            "italic(d) >= 0", 
                            "italic(d) >= 0 ~~~ italic(p) >= 0.25")),
  core_size = factor(c, labels = paste("italic(c) ==", cvals, sep = '')),
  acc = factor(acc, levels = c("mean_fp", "mean_tn", "mean_tp", "mean_fn"),
               labels = c("False core", "True periphery", "True core", "False periphery"))
)]

# Function to create ggplot for specified models
create_plot <- function(data, models, facet_by = c("core_size", "type")) {
  # Define facet formula based on the selected facets
  facet_formula <- if (length(facet_by) == 1) {
    as.formula(paste("model ~", facet_by[1]))
  } else {
    as.formula(paste("model ~", paste(facet_by, collapse = " + ")))
  }
  
  # Plot creation
  ggplot(data = data[model %in% models], aes(x = k, y = value, fill = acc)) +
    geom_bar(stat = "identity") +
    geom_hline(aes(yintercept = c)) +
    geom_text(aes(x = 2, y = c + 2.5, label = "periphery"), size = 2.5) +
    geom_text(aes(x = 2, y = c - 2.5, label = "core"), size = 2.5) +
    facet_grid(facet_formula, labeller = label_parsed) +
    labs(y = "Average confusion matrix", color = "", fill = "", linetype = "") +
    scale_fill_manual(values = c("red", "green2", "chartreuse", "red3")) +
    scale_x_continuous(breaks = seq(1, 3, by = 1)) +
    theme(legend.position = "top",
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank())
}

# Create plots
p1.1 <- create_plot(data_long, c("italic(d) == 'NA'", "italic(d) >= 0"),
                    facet_by = "core_size")
p1.2 <- create_plot(data_long, c("italic(d) == 'NA'", "italic(d) >= 0 ~~~ italic(p) >= 0.25"),
                    facet_by = "core_size")

# Function to prepare data for modeling and plotting
prepare_model_data <- function(data, model_names) {
  data <- data[name %in% model_names]
  data[, name := factor(name, levels = model_names)]
  data[, ntw := rep(1:(nrow(data) / 2), each = 2)]
  data[, ntw := as.factor(ntw)]
  data[, c := as.factor(c)]
  data[, success := TP + TN]
  data[, fail := FP + FN]
  return(data)
}

# Model and plot function ----
run_model_and_plot <- function(data, model_name) {
  model <- glmmTMB(cbind(success, fail) ~ name * k * c + (1 | ntw),
                   family = binomial(link = 'logit'), data = data)
  predvals <- ggpredict(model, terms = c("k [all]", "name", "c"),
                        ci_level = .995) # 99.5% Confidence intervals
  predvals <- data.table(predvals)
  
  ggplot(data = predvals, aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high,
                              color = group, fill = group)) +
    geom_ribbon(alpha = .5, linewidth = .1) + 
    geom_line(linewidth = .25) +
    scale_color_manual(values = c("red3", "navyblue"),
                       labels = c(expression(italic(d) == "NA"), expression(italic(d) >= 0))) +
    scale_fill_manual(values = c("orange", "royalblue"),
                      labels = c(expression(italic(d) == "NA"), expression(italic(d) >= 0))) +
    facet_grid("Confidence intervals (99.5%)" ~ facet) +
    labs(x = expression(italic(k)), y = "Actor assignment accuracy", color = "", fill = "") +
    scale_y_continuous(labels = percent_format(scale = 100)) +
    scale_x_continuous(breaks = seq(1, 3, by = 1)) +
    theme(strip.text.x = element_blank(), legend.position = "top")
}

# Run analysis and plotting for each model
data1 <- prepare_model_data(data, c("ucinet", "minden"))
p2.1 <- run_model_and_plot(data1, "minden")

# Combined plots
combined_plot1 <- p1.1 + p2.1 + plot_layout(ncol = 1, heights = c(2, 1))

# Save the combined plot
tiff(filename = "Fig7.1.tiff", width = 27, height = 20, 
     units = "cm", compression = "lzw", bg = "white", res = 1000)
combined_plot1
dev.off()

# Run analysis and plotting for the second model
data2 <- prepare_model_data(data, c("ucinet", "pcore"))
p2.2 <- run_model_and_plot(data2, "pcore")

# Combined plots for the second analysis
combined_plot2 <- p1.2 + p2.2 + plot_layout(ncol = 1, heights = c(2, 1))

# Save the combined plot
tiff(filename = "Fig7.2.tiff", width = 27, height = 20,
     units = "cm", compression = "lzw", bg = "white", res = 1000)
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
  lapply(simntw[sample],cp.pcore,delta=0,p=0.25),
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
results3 <- lapply(simntw,cp.pcore,delta=0,p=0.25) # min-density blocks, plus p-core 0.25

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

# Data wrangling for visualization ----
means <- data[, lapply(.SD, mean), .SDcols = c("TN", "FP", "FN", "TP"), by = .(name, k, c, conn)]
setnames(means, c("TN", "FP", "FN", "TP"), paste0("mean_", c("tn", "fp", "fn", "tp")))

# Reshape data to long format
data_long <- as.data.table(pivot_longer(as_tibble(means), cols = starts_with("mean_"), names_to = "acc", values_to = "value"))

# Rename factors
data_long[, `:=`(
  model = factor(name, levels = c("ucinet", "minden", "pcore"),
                 labels = c("italic(d) == 'NA'", 
                            "italic(d) >= 0", 
                            "italic(d) >= 0 ~~~ italic(p) >= 0.25")),
  core_size = factor(c, labels = paste("italic(c) ==", cvals, sep = '')),
  type = factor(conn,levels=c('no','core','periphery'),
                labels=c('Unconnected','Core','Periphery')),
  acc = factor(acc, levels = c("mean_fp", "mean_tn", "mean_tp", "mean_fn"),
               labels = c("False core", "True periphery", "True core", "False periphery"))
)]

# Create plots
p3.1 <- create_plot(data_long, c("italic(d) == 'NA'", "italic(d) >= 0"),
                    facet_by = c("core_size", "type"))
p3.2 <- create_plot(data_long, c("italic(d) == 'NA'", "italic(d) >= 0 ~~~ italic(p) >= 0.25"),
                    facet_by = c("core_size", "type"))

# Function to prepare data for modeling and plotting
prepare_model_data <- function(data, model_names) {
  data <- data[name %in% model_names]
  data[, name := factor(name, levels = model_names)]
  data[, ntw := rep(1:(nrow(data) / 2), each = 2)]
  data[, ntw := as.factor(ntw)]
  data[, c := as.factor(c)]
  data[, conn := factor(conn,levels=c('no','core','periphery'))] 
  data[, success := TP + TN]
  data[, fail := FP + FN]
  return(data)
}

# Model and plot function ----
run_model_and_plot <- function(data, model_name) {
  model <- glmmTMB(cbind(success, fail) ~ name * k * c * conn + (1 | ntw),
                   family = binomial(link = 'logit'), data = data)
  predvals <- ggpredict(model, terms = c("k [all]", "name", "conn", "c"),
                        ci_level = .995) # 99.5% Confidence intervals
  predvals <- data.table(predvals)
  
  ggplot(data = predvals, aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high,
                              color = group, fill = group)) +
    geom_ribbon(alpha = .5, linewidth = .1) + 
    geom_line(linewidth = .25) +
    scale_color_manual(values = c("red3", "navyblue"),
                       labels = c(expression(italic(d) == "NA"), expression(italic(d) >= 0))) +
    scale_fill_manual(values = c("orange", "royalblue"),
                      labels = c(expression(italic(d) == "NA"), expression(italic(d) >= 0))) +
    facet_grid("Confidence intervals (99.5%)" ~ panel + facet) +
    labs(x = expression(italic(k)), y = "Actor assignment accuracy", color = "", fill = "") +
    scale_y_continuous(labels = percent_format(scale = 100)) +
    scale_x_continuous(breaks = seq(1, 3, by = 1)) +
    theme(strip.text.x = element_blank(), legend.position = "top")
}

# Run analysis and plotting for each model
data1 <- prepare_model_data(data, c("ucinet", "minden"))
p4.1 <- run_model_and_plot(data1, "minden")

# Combined plots
combined_plot1 <- p3.1 + p4.1 + plot_layout(ncol = 1, heights = c(2, 1))

# Save the combined plot
tiff(filename = "Fig9.1.tiff", width = 27, height = 20, 
     units = "cm", compression = "lzw", bg = "white", res = 1000)
combined_plot1
dev.off()

# Run analysis and plotting for the second model
data2 <- prepare_model_data(data, c("ucinet", "pcore"))
p4.2 <- run_model_and_plot(data2, "pcore")

# Combined plots for the second analysis
combined_plot2 <- p3.2 + p4.2 + plot_layout(ncol = 1, heights = c(2, 1))

# Save the combined plot
tiff(filename = "Fig7.2.tiff", width = 27, height = 20,
     units = "cm", compression = "lzw", bg = "white", res = 1000)
combined_plot2
dev.off()

################################################################################