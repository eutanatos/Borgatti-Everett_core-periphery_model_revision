########################################################################################################################
## REVISING THE BORGATTI-EVERETT CORE-PERIPHERY MODEL
## (1) Illustration figures
## R script written by José Luis Estévez (University of Helsinki / Vaestoliitto)
## Date: Aug 3rd, 2024
########################################################################################################################

# R PACKAGES REQUIRED
library(data.table);library(ggplot2);library(ggpubr)

# CLEAN ENVIRONMENT
rm(list=ls())

# Set theme for plots
theme_set(theme_minimal() + theme(legend.position = 'none'))

########################################################################################################################

# Let's create a network for visualization purposes
data <- data.table(Var1 = rep(1:10,time=10),
                   Var2 = rep(1:10,each=10),
                   value = c(NA,1,1,0,1,0,1,1,0,1,
                             0,NA,0,1,1,0,0,1,1,0,
                             1,1,NA,1,1,1,0,0,0,0,
                             1,0,0,NA,0,0,1,0,0,0,
                             0,1,0,1,NA,0,1,0,0,0,
                             0,0,0,0,0,NA,0,0,0,0,
                             0,0,0,0,1,0,NA,1,0,1,
                             0,1,0,0,1,0,1,NA,0,0,
                             0,0,0,1,0,0,0,0,NA,0,
                             1,1,0,0,0,0,1,1,0,NA))

data[,Var1 := factor(Var1)]
data[,Var2 := factor(Var2,levels=rev(1:10))]

p1 <- ggplot(data = data, aes(x = Var1, y = Var2)) +
  geom_tile(fill = ifelse(data$Var1 == data$Var2,'grey40','white')) +
  # borders of the matrix
  geom_vline(xintercept = c(0.5,10.5), color = 'black') +
  geom_hline(yintercept = c(0.5,10.5), color = 'black') +
  # border of the cells
  geom_vline(xintercept = seq(1.5,9.5,by=1), color = 'grey80') +
  geom_hline(yintercept = seq(1.5,9.5,by=1), color = 'grey80') +
  # partition of the blockmodel
  geom_vline(xintercept = 5.5, color = 'red',linetype='dashed',linewidth=1) +
  geom_hline(yintercept = 5.5, color = 'red',linetype='dashed',linewidth=1) +
  labs(x = '', y = '') +
  scale_x_discrete(position = "top") +
  geom_text(aes(label = ifelse(data$value == 1,value,''))) +
  annotate("text",x=3,y=8,label=expression(bold(A)^cc),color='navyblue',size=15,alpha=0.25,family='serif') +
  annotate("text",x=8,y=8,label=expression(bold(A)^cp),color='navyblue',size=15,alpha=0.25,family='serif') +
  annotate("text",x=8,y=3,label=expression(bold(A)^pp),color='navyblue',size=15,alpha=0.25,family='serif') +
  annotate("text",x=3,y=3,label=expression(bold(A)^pc),color='navyblue',size=15,alpha=0.25,family='serif') 

# Ideal pattern with delta
data[,ip1 := paste('italic(d)')]
core <- 1:5
periphery <- 6:10
data[Var1 %in% core & Var2 %in% core,ip1 := 1]
data[Var1 %in% periphery & Var2 %in% periphery,ip1 := 0]
data[Var1 == Var2,ip1 := NA] # remove diagonal

p2 <- ggplot(data = data, aes(x = Var1, y = Var2)) +
  geom_tile(fill = ifelse(data$Var1 == data$Var2,'grey40','white')) +
  # borders of the matrix
  geom_vline(xintercept = c(0.5,10.5), color = 'black') +
  geom_hline(yintercept = c(0.5,10.5), color = 'black') +
  # border of the cells
  geom_vline(xintercept = seq(1.5,9.5,by=1), color = 'grey80') +
  geom_hline(yintercept = seq(1.5,9.5,by=1), color = 'grey80') +
  # partition of the blockmodel
  geom_vline(xintercept = 5.5, color = 'red',linetype='dashed',linewidth=1) +
  geom_hline(yintercept = 5.5, color = 'red',linetype='dashed',linewidth=1) +
  labs(x = '', y = '') +
  scale_x_discrete(position = "top") +
  geom_text(aes(label = ifelse(data$ip1 != 0,ip1,'')),parse=TRUE) +
  annotate("text",x=3,y=8,label=expression(bold(B)(bold(v))^cc),color='navyblue',size=15,alpha=0.25,family='serif') +
  annotate("text",x=8,y=8,label=expression(bold(B)(bold(v))^cp),color='navyblue',size=15,alpha=0.25,family='serif') +
  annotate("text",x=8,y=3,label=expression(bold(B)(bold(v))^pp),color='navyblue',size=15,alpha=0.25,family='serif') +
  annotate("text",x=3,y=3,label=expression(bold(B)(bold(v))^pc),color='navyblue',size=15,alpha=0.25,family='serif') 

# Visualization
tiff(filename="Fig1.tiff",
     width=22, height=11,units="cm", 
     compression="lzw",bg="white",res=1000
)
ggarrange(p1,p2,ncol=2,labels=c('A','B'))
dev.off()

########################################################################################################################

# P-CORE ROUTINE

# Single out the core
core <- data[Var1 %in% 1:5 & Var2 %in% 1:5]
core <- core[order(rev(Var2),Var1)]

c1 <- ggplot(data = core,
             aes(x = Var1, y = Var2)) +
  geom_tile(fill=ifelse(core$Var1 == core$Var2,'grey40','white')) +
  # borders of the matrix
  geom_vline(xintercept = c(0.5,5.5), color = 'black') +
  geom_hline(yintercept = c(0.5,5.5), color = 'black') +
  # border of the cells
  geom_vline(xintercept = seq(1.5,4.5,by=1), color = 'grey80') +
  geom_hline(yintercept = seq(1.5,4.5,by=1), color = 'grey80') +
  labs(x = '', y = '') +
  scale_x_discrete(position = "top") +
  geom_text(aes(label = ifelse(core$value == 1,1,''))) +
  annotate("text",x=3,y=3,label=expression(bold(A)^cc),
           color='navyblue',size=15,alpha=0.25,family='serif')

# Rearranged the values
core[,value2 := c(1,1,1,0,NA,
                  1,1,0,0,NA,
                  1,1,1,1,NA,
                  1,0,0,0,NA,
                  1,1,0,0,NA)]

c2 <- ggplot(data = core,
             aes(x = Var1, y = Var2)) +
  geom_tile(fill=ifelse(is.na(core$value2),'grey40','white')) +
  # borders of the matrix
  geom_vline(xintercept = c(0.5,5.5), color = 'black') +
  geom_hline(yintercept = c(0.5,5.5), color = 'black') +
  # border of the cells
  geom_vline(xintercept = seq(1.5,4.5,by=1), color = 'grey80') +
  geom_hline(yintercept = seq(1.5,4.5,by=1), color = 'grey80') +
  labs(x = '', y = '') +
  scale_x_discrete(position = "top") +
  geom_text(aes(label = ifelse(core$value2 == 1,1,''))) +
  annotate("text", x=3, y=3,
           label=expression(bold(A)^{cc~"\u2192"}),
           color='navyblue',size=15,alpha=0.25,family='serif') 

core[,value3 := c(1,1,1,1,1,
                  1,1,0,1,1,
                  0,1,0,1,1,
                  0,0,0,0,0,
                  NA,NA,NA,NA,NA)]

c3 <- ggplot(data = core,
             aes(x = Var1, y = Var2)) +
  geom_tile(fill=ifelse(is.na(core$value3),'grey40','white')) +
  # borders of the matrix
  geom_vline(xintercept = c(0.5,5.5), color = 'black') +
  geom_hline(yintercept = c(0.5,5.5), color = 'black') +
  # border of the cells
  geom_vline(xintercept = seq(1.5,4.5,by=1), color = 'grey80') +
  geom_hline(yintercept = seq(1.5,4.5,by=1), color = 'grey80') +
  labs(x = '', y = '') +
  scale_x_discrete(position = "top") +
  geom_text(aes(label = ifelse(core$value3 == 1,1,''))) +
  annotate("text", x=3, y=3,
           label=expression(bold(A)^{cc~"\u2193"}),
           color='navyblue',size=15,alpha=0.25,family='serif')

# Ideal blocks
core[,ideal1 := "1"]

for(i in 1:nrow(core)){
  if(core$Var1[i] %in% 3:4){
    core$ideal1[i] <- paste('italic(a)[',core$Var2[i],core$Var1[i],']^',
                         expression(paste("cc\u2192")),sep='')
  }
}

core[Var1 == 5]$ideal1 <- NA # add the diagonal at the end

c4 <- ggplot(data = core,
             aes(x = Var1, y = Var2)) +
  geom_tile(fill=ifelse(is.na(core$ideal1),'grey40','white')) +
  # borders of the matrix
  geom_vline(xintercept = c(0.5,5.5), color = 'black') +
  geom_hline(yintercept = c(0.5,5.5), color = 'black') +
  # border of the cells
  geom_vline(xintercept = seq(1.5,4.5,by=1), color = 'grey80') +
  geom_hline(yintercept = seq(1.5,4.5,by=1), color = 'grey80') +
  # Cut-off point (0.5 p-core)
  geom_vline(xintercept = 2.5, color = 'orange3',linetype='dashed',linewidth=1) +
  labs(x = '', y = '') +
  scale_x_discrete(position = "top") +
  geom_text(aes(label = ideal1),parse=TRUE) +
  annotate("text", x=3, y=3,
           label=expression(bold(B)(bold(v))^{cc~"\u2192"}),
           color='navyblue',size=15,alpha=0.25,family='serif')

# second version
core[,ideal2 := "1"]

for(i in 1:nrow(core)){
  if(core$Var2[i] %in% 3:4){
    core$ideal2[i] <- paste('italic(a)[',core$Var2[i],core$Var1[i],']^',
                            expression(paste("cc\u2193")),sep='')
  }
}

core[Var2 == 5]$ideal2 <- NA # add the diagonal at the end

c5 <- ggplot(data = core,
       aes(x = Var1, y = Var2)) +
  geom_tile(fill=ifelse(is.na(core$ideal2),'grey40','white')) +
  # borders of the matrix
  geom_vline(xintercept = c(0.5,5.5), color = 'black') +
  geom_hline(yintercept = c(0.5,5.5), color = 'black') +
  # border of the cells
  geom_vline(xintercept = seq(1.5,4.5,by=1), color = 'grey80') +
  geom_hline(yintercept = seq(1.5,4.5,by=1), color = 'grey80') +
  # Cut-off point (0.5 p-core)
  geom_hline(yintercept = 3.5, color = 'orange3',linetype='dashed',linewidth=1) +
  labs(x = '', y = '') +
  scale_x_discrete(position = "top") +
  geom_text(aes(label = ideal2),parse=TRUE) +
  annotate("text", x=3, y=3,
           label=expression(bold(B)(bold(v))^{cc~"\u2193"}),
           color='navyblue',size=15,alpha=0.25,family='serif')

# Visualization
empty_plot <- ggplot() + theme_void()

tiff(filename="Fig2.tiff",
     width=24, height=16,units="cm", 
     compression="lzw",bg="white",res=1000
)
ggarrange(ggarrange(empty_plot,c1,empty_plot,ncol=1,labels=c('','A',''),heights = c(1,2,1)),
          ggarrange(c2,c4,c3,c5,nrow=2,ncol=2,labels=c('B','C','D','E')),
          ncol=2,widths = c(1,2))
dev.off()

########################################################################################################################