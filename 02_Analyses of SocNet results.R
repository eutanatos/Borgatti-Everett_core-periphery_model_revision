################################################################################
## REVISING THE BORGATTI-EVERETT CORE-PERIPHERY MODEL
## (2) Analyses of SocNet results
## R script written by José Luis Estévez (University of Helsinki / Vaestoliitto)
## Date: Aug 4th, 2024
################################################################################

# R PACKAGES REQUIRED
library(data.table);library(ggplot2);library(purrr);library(stringr);library(ggpubr);library(sna);library(igraph)

# CLEAN ENVIRONMENT
rm(list=ls())

# Set theme for plots
theme_set(theme_bw())

################################################################################

# DATA LOADING ----

# befig1, Baker, Galtung's and Zachary's networks
mtx <- list()
mtx$befig1 <- read.table('networks/befig1.txt',header=TRUE,sep='\t')
mtx$baker <- read.table('networks/baker.txt',header=TRUE,sep='\t')
mtx$galtung <- read.table('networks/galtung.txt',header=TRUE,sep='\t')
mtx$zachary <- read.table('networks/zachary.txt',header=TRUE,sep='\t')

# Remember to remove the first column (the column names), and assign column and row names
for(i in 1:length(mtx)){
  colnames(mtx[[i]])[-1] <- rownames(mtx[[i]]) <- toupper(mtx[[i]][,1]) # upper cases
  mtx[[i]] <- as.matrix(mtx[[i]][,-1]) # remove first column (which are row names)
  diag(mtx[[i]]) <- NA # remove the diagonal
}

# Socnet.se output
socnet <- data.table(read.table('socnet output.csv',header=TRUE,sep=','))
# core size
socnet[,coresize := str_count(socnet$partition,'0')]

################################################################################

# NETWORK VISUALIZATIONS ----
grphs <- lapply(mtx,graph_from_adjacency_matrix,mode='undirected') # turn into igraph objects

# layouts
set.seed(0708)
ltys <- lapply(grphs,layout_with_kk)

# Network names
ntwnames <- c('BEfig1','Baker','Galtung','Zachary')

# Visualization
par(mfrow=c(2,2))
for(i in 1:4){
  plot.igraph(grphs[[i]],
              vertex.color='grey95',vertex.size=23,vertex.frame.color='black',
              vertex.label.family = 'Arial',vertex.label.cex=.75,vertex.label.font=2,vertex.label.color='black',
              edge.color='grey70',edge.width=2.5,
              layout = ltys[[i]],
              main=ntwnames[i])
}

#############################

# Add core members (based on standard Borgatti-Everett method)
stpart <- socnet[pcore=='com' & intercategorical=='dnc',.(network,partition)]
# stpart <- socnet[pcore=='com' & intercategorical=='denuci' & d %in% c(0.3333,0.3846,0.25,0.2044),.(network,partition)]
#stpart <- socnet[pcore=='com' & intercategorical=='denmin' & d == 0,.(network,partition)]
coremembers <- list() 

for(i in 1:nrow(stpart)){
  ntwname <- stpart[i,]$network
  pttn <- stpart[i,]$partition # take the partition
  pttn <- gsub('(','',pttn,fixed=TRUE) # remove the parentheses
  pttn <- gsub(')','',pttn,fixed=TRUE) 
  pttn <- strsplit(pttn,split=';') # turn into vector
  core <- rownames(mtx[[i]])[pttn[[1]] == 0] # members of the core
  coremembers[[paste(ntwname)]] <- core
}

# Find periphery-periphery ties, if there are any
trnsf <- purrr::compose(as.data.table,as_edgelist) # turn into edgelist/data.table
edgelists <- lapply(grphs,trnsf)

for(i in 1:length(edgelists)){
  E(grphs[[i]])$periperi <- edgelists[[i]][,!(V1 %in% coremembers[[i]]) & 
                                             !(V2 %in% coremembers[[i]])]
}

# Goodness of fit values
gof <- socnet[pcore=='com' & intercategorical=='dnc',.(network,gof)]

# Visualization
tiff(filename="Fig3.tiff",
     width=25, height=25,units="cm", 
     compression="lzw",bg="white",res=1000
)
par(mfrow=c(2,2))
for(i in 1:4){
  plot.igraph(grphs[[i]],
              vertex.color=ifelse(V(grphs[[i]])$name %in% coremembers[[i]],'grey10','grey95'),
              vertex.size=23,vertex.frame.color='black',
              vertex.label.family = 'Arial',vertex.label.cex=.75,vertex.label.font=2,
              vertex.label.color=ifelse(V(grphs[[i]])$name %in% coremembers[[i]],'white','black'),
              edge.color='grey70',
              edge.width=2.5,
              layout = ltys[[i]],
              mark.groups = coremembers[i],mark.col='grey90',mark.border=NA,
              main=ntwnames[i])
  # add GOF
  text(x=-1.1, y=-1.1, labels=bquote(paste(italic(r),' = ',.(round(gof$gof[i],4)),sep='')), cex=1.25,col='red3')
}
dev.off()

################################################################################

# SOME NETWORK DESCRIPTIVES ----

# Combined function to obtained average degree
# Define the functions
func1 <- function(x) rowSums(x,na.rm=TRUE)
func2 <- function(x) mean(x)
# Combine the functions
combined_func <- function(x) {
  result1 <- func1(x)
  result2 <- func2(result1)
  return(result2)
}

desc <- data.table(network = names(grphs),
           nodes = unlist(lapply(grphs,vcount)),
           ties = unlist(lapply(grphs,ecount)),
           density = unlist(lapply(mtx,gden)),
           ave.degree = unlist(lapply(mtx,combined_func)),
           transitivity = unlist(lapply(mtx,gtrans,measure='weak')),
           diameter = unlist(lapply(grphs,diameter)))

desc
write.table(desc,'descriptives.csv',sep=',',row.names = FALSE)

################################################################################

# INTER-CATEGORICAL DENSITIES (OUTPUT) ----

# Let's see inter-categorical blocks' densities
# Since these are all undirected networks, we can just use one block
interblock <- mtx
for(i in 1:length(interblock)){
  members <- row.names(mtx[[i]]) %in% coremembers[[i]]
  interblock[[i]] <- interblock[[i]][members,!members]
}

# Inter-categorical densities
intercat <- data.table(ties = unlist(lapply(interblock,sum)), # ties
                       cells = unlist(lapply(interblock,length))) # possible ties
dobserved <- data.table(network = names(grphs),
                        doutput = round(intercat[,ties/cells],4)) # values used in socnet.se
dobserved

################################################################################

# INTER-CATEGORICAL DENSITY SEARCH ----
icds <- socnet[intercategorical %in% c('denuci','den','denmin') 
               & pcore == 'com'
               & as.character(d) %in% seq(0,1,by=0.05)]
icds[,network := factor(network,levels=c('befig1','baker','galtung','zachary'),
                        labels=c('BEfig1','Baker','Galtung','Zachary'))]
icds[,intercategorical := factor(intercategorical,levels = c('denuci','den','denmin'),
                                 labels=c("Standard method","Exact-density","Minimum-density"))]

dobserved[,network := factor(network,levels=c('befig1','baker','galtung','zachary'),
                        labels=c('BEfig1','Baker','Galtung','Zachary'))]

# Visualization
tiff(filename="Fig4.tiff",
     width=20, height=15,units="cm", 
     compression="lzw",bg="white",res=1000
)
ggplot(data=icds[solution == 1],aes(x=d,y=gof)) +
  #geom_bar(stat='identity',fill='grey90',color='grey50') + 
  geom_line(lty='solid',alpha=.5) + geom_point() + 
  geom_vline(data=dobserved,aes(xintercept = doutput),color='red3',lty='dashed') +
  geom_text(data=dobserved,aes(x=doutput-0.15,y=.2,label=doutput),color='red3',angle=90,vjust=1.5,size=3) +
  ylim(c(0.1,1)) +
  labs(x=expression(paste(italic(d),' (inter-categorical density)',sep='')),
       y=expression(paste(italic('r'),' (goodness-of-fit measure)',sep=''))) +
  facet_grid(intercategorical~network) 
dev.off()

################################################################################

# VISUALIZATION OF DISTORTING EFFECT ON BEFIG1 ----

be <- as.data.table(reshape2::melt(mtx$befig1))

# To identify the block of each pair.
be[,block := ifelse(Var1 %in% 1:4 & Var2 %in% 1:4,'Core',
                    ifelse(Var1 %in% 5:10 & Var2 %in% 5:10,'Periphery','Inter-categorical'))]

# Ideal pattern matrix (standard BE method)
be[,ideal1 := ifelse(block == 'Core',1,ifelse(block == 'Periphery',0,1/3))] # d = 1/3
be[,ideal2 := ifelse(ideal1 == 1/3,0.2833,ideal1)] # d = 0.2833
be[,ideal3 := ifelse(ideal1 == 1/3,0.3833,ideal1)] # d = 0.3833

# Ideal pattern matrix (using exact-density ideal blocks)
be[,ideal4 := value] # d = 0.2833
be[Var1 == 1 & Var2 == 5]$ideal4 <- be[Var1 == 5 & Var2 == 1]$ideal4 <- 0

be[,ideal5 := value] # d = 0.3833
be[Var1 == 1 & Var2 == 6]$ideal5 <- be[Var1 == 6 & Var2 == 1]$ideal5 <- 1

# Correlation tests
# Standard method
round(cor.test(be$value,be$ideal1)$estimate,4)
round(cor.test(be$value,be$ideal2)$estimate,4)
round(cor.test(be$value,be$ideal3)$estimate,4)
# Exact density
round(cor.test(be$value,be$value)$estimate,4)
round(cor.test(be$value,be$ideal4)$estimate,4)
round(cor.test(be$value,be$ideal5)$estimate,4)

# Visualization 
be33 <- ggplot(data=be,aes(x=value,y=ideal1)) +
  geom_smooth(method='lm',se=FALSE,color='grey20') +
  geom_jitter(aes(fill=block,shape=block),width = 0.025, height = 0.025,alpha=2/3) +
  annotate("text",x=.2,y=.9,label = parse(text = "italic(d) == 0.3333")) +
  annotate("text",x=.2,y=.8,label = parse(text = "italic(r) == 0.6686")) +  
  scale_shape_manual(values = c('Core'=24,'Inter-categorical'=23,'Periphery'=25)) +
  labs(x='',y='Ideal',fill='Block',shape='Block') 
be28 <- ggplot(data=be,aes(x=value,y=ideal2)) +
  geom_smooth(method='lm',se=FALSE,color='grey20') +
  geom_jitter(aes(fill=block,shape=block),width = 0.025, height = 0.025,alpha=2/3) +
  annotate("text",x=.2,y=.9,label = parse(text = "italic(d) == 0.2833")) +
  annotate("text",x=.2,y=.8,label = parse(text = "italic(r) == 0.6664")) +  
  scale_shape_manual(values = c('Core'=24,'Inter-categorical'=23,'Periphery'=25)) +
  labs(x='',y='Ideal',fill='Block',shape='Block')
be38 <- ggplot(data=be,aes(x=value,y=ideal3)) +
  geom_smooth(method='lm',se=FALSE,color='grey20') +
  geom_jitter(aes(fill=block,shape=block),width = 0.025, height = 0.025,alpha=2/3) +
  annotate("text",x=.2,y=.9,label = parse(text = "italic(d) == 0.3833")) +
  annotate("text",x=.2,y=.8,label = parse(text = "italic(r) == 0.6665")) +  
  scale_shape_manual(values = c('Core'=24,'Inter-categorical'=23,'Periphery'=25)) +
  labs(x='Observed',y='Ideal',fill='Block',shape='Block') 

edb1 <- ggplot(data=be,aes(x=value,y=value)) +
  geom_smooth(method='lm',se=FALSE,color='grey20') +
  geom_jitter(aes(fill=block,shape=block),width = 0.025, height = 0.025,alpha=2/3) +
  annotate("text",x=.2,y=.9,label = parse(text = "italic(d) == 0.3333")) +
  annotate("text",x=.2,y=.8,label = parse(text = "italic(r) == 1")) +  
  scale_shape_manual(values = c('Core'=24,'Inter-categorical'=23,'Periphery'=25)) +
  labs(x='',y='',fill='Block',shape='Block') 
edb2 <- ggplot(data=be,aes(x=value,y=ideal4)) +
  geom_smooth(method='lm',se=FALSE,color='grey20') +
  geom_jitter(aes(fill=block,shape=block),width = 0.025, height = 0.025,alpha=2/3) +
  annotate("text",x=.2,y=.9,label = parse(text = "italic(d) == 0.2833")) +
  annotate("text",x=.2,y=.8,label = parse(text = "italic(r) == 0.9484")) +  
  scale_shape_manual(values = c('Core'=24,'Inter-categorical'=23,'Periphery'=25)) +
  labs(x='',y='',fill='Block',shape='Block') 
edb3 <- ggplot(data=be,aes(x=value,y=ideal5)) +
  geom_smooth(method='lm',se=FALSE,color='grey20') +
  geom_jitter(aes(fill=block,shape=block),width = 0.025, height = 0.025,alpha=2/3) +
  annotate("text",x=.2,y=.9,label = parse(text = "italic(d) == 0.3833")) +
  annotate("text",x=.2,y=.8,label = parse(text = "italic(r) == 0.9504")) +  
  scale_shape_manual(values = c('Core'=24,'Inter-categorical'=23,'Periphery'=25)) +
  labs(x='Observed',y='',fill='Block',shape='Block')

# Visualization
tiff(filename="Fig5.tiff",
     width=14, height=21,units="cm", 
     compression="lzw",bg="white",res=1000
)
ggarrange(be33,edb1,be28,edb2,be38,edb3,
          nrow=3,ncol=2,labels=c('A','B','C','D','E','F'),common.legend = TRUE)
dev.off()

################################################################################

# P-CORE ----

# Let's see the effect of a smaller p-core value
forplot <- socnet[intercategorical %in% c('dnc','denmin')] 
forplot <- forplot[is.na(d) | (!is.na(d) & d == 0)]

forplot[pcore == 'com',pcore := 1] # A complete core equals 1
forplot[,network := factor(network,levels=c('befig1','baker','galtung','zachary'),
                           labels=c('BEfig1','Baker','Galtung','Zachary'))]
forplot[,intercategorical := factor(intercategorical,
                                    levels=c('dnc','denmin'),
                                    labels=c('italic(d)=="NA"','italic(d)>=0'))]
forplot[,pcore := as.numeric(pcore)]
# Number of solutions
forplot <- forplot[,length(solution),by=.(network,size,intercategorical,pcore,gof,coresize)]

# Densities of the full networks (to show in figure)
ntwdens <- data.table(network = names(grphs),
                      density = unlist(lapply(grphs,edge_density)))
ntwdens[,network := factor(network,levels=c('befig1','baker','galtung','zachary'),
                           labels=c('BEfig1','Baker','Galtung','Zachary'))]

# Visualization
tiff(filename="Fig6.tiff",
     width=20, height=15,units="cm", 
     compression="lzw",bg="white",res=1000
)
ggplot(data=forplot) +
  geom_vline(data=ntwdens,aes(xintercept=density),linetype='dashed',color='orange3') +
  geom_rect(data=ntwdens,aes(xmin=density,xmax=0,ymin=0,ymax=1),fill='bisque1',color=NA,alpha=.5) +
  geom_text(data=ntwdens,aes(x=density,y=.5,label=paste("Network density:",round(density,3),sep=' ')),color='orange3',angle=90,vjust=1.5,size=3) +
  geom_point(aes(x=pcore,y=coresize/size),shape=21,color='black',fill='navyblue',size=3.5) +
  geom_text(aes(x=pcore,y=coresize/size,group=network,label=V1),color='white',size=2) +
  geom_line(aes(x=pcore,y=gof,group=network),color='red3',alpha=.5) +
  geom_point(aes(x=pcore,y=gof),shape=23,color='black',fill='red3',size=2.5) +
  scale_y_continuous(name=parse(text="italic(r[w])"),
                     sec.axis = sec_axis(~./1,name="Core members (proportion)")) +
  theme(axis.title.y.left=element_text(color="red3"),
        axis.title.y.right=element_text(color='navyblue')) +
  scale_x_reverse() + labs(x=expression(paste(italic(p),'-core value',sep=''))) +
  facet_grid(intercategorical~network,labeller=label_parsed) 
dev.off()

################################################################################

# Let's extend the core periphery solutions to add core members with p-cores values >= 0.5
# Alternative cores 
altcores <- socnet[intercategorical == 'denmin' 
                   & pcore %in% c(0.5,0.6,0.7,0.8,0.9,'com')]
altcores <- altcores[is.na(d) | (!is.na(d) & d == 0)]

altcoremembers <- list() 

for(i in 1:nrow(altcores)){
  ntwname <- altcores[i,]$network
  pttn <- altcores[i,]$partition # take the partition
  pttn <- gsub('(','',pttn,fixed=TRUE) # remove the parentheses
  pttn <- gsub(')','',pttn,fixed=TRUE) 
  pttn <- strsplit(pttn,split=';') # turn into vector
  core <- rownames(mtx[[ntwname]])[pttn[[1]] == 0] # members of the core
  altcoremembers[[paste(ntwname,i)]] <- core
}

# Let's divide the information by network
pseudocores <- vector('list',4)
names(pseudocores) <- tolower(ntwnames)

for(i in 1:length(pseudocores)){
  ntwname <- names(pseudocores)[i]
  pseudocores[[i]] <- altcoremembers[startsWith(names(altcoremembers),ntwname)]
}

# Now we can obtain the unions and intersections
pseudocoresunion <- pseudocoresintersect <- pseudocores
for(i in 1:length(pseudocores)){
  pseudocoresunion[[i]] <- Reduce(union,pseudocores[[i]])
  pseudocoresintersect[[i]] <- Reduce(intersect,pseudocores[[i]])
}

# Visualization
par(mfrow=c(2,2))
for(i in 1:4){
  plot.igraph(grphs[[i]],
              vertex.color=ifelse(V(grphs[[i]])$name %in% pseudocoresintersect[[i]],'grey10',
                                  ifelse(V(grphs[[i]])$name %in% pseudocoresunion[[i]],'grey50','grey95')),
              vertex.size=23,vertex.frame.color='black',
              vertex.label.family = 'Arial',vertex.label.cex=.75,vertex.label.font=2,
              vertex.label.color=ifelse(V(grphs[[i]])$name %in% pseudocoresunion[[i]],'white','black'),
              edge.color='grey70',
              edge.width=2.5,
              layout = ltys[[i]],
              mark.groups = list(pseudocoresunion[[i]],pseudocoresintersect[[i]]),
              mark.col=c('gray95','grey85'),
              mark.border = c('grey70',NA),
              main=ntwnames[i])
  # add legend
  text(x=-1.1, y=-1.1, labels=bquote(italic(p) >= 0.5), cex=1.25,col='red3')
}

################################################################################