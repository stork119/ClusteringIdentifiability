#########################################
### Sesnsitivity Analysis Clustering  ###
#########################################

### Arguments :
# File directory
# Sensitivity Matrix File Directory
# Output Directories 
# zeta
# delta
# Parameters Names File Directory/Directories 
###

# Clear workspace
rm(list=ls()) 
# LIBRARIES
source("SACluster.R")

#### arguments####

delta            <- 0.95
zeta             <- 1

directory.models <- '../../Models/'

#### NFKB EXAMPLE ####
directory.folder <- paste(directory.models,
                          'NFkB/Output/2015-07-17/t-175-tdeg-0/',
                          sep = "")
directory.SM     <- paste(directory.folder, 'sense.csv', sep = "")
directory.FIM     <- paste(directory.folder, 'FIM.csv', sep = "")
directory.output <- paste(directory.folder,
                          "clustering/",
                          'delta/',
                          '/', 
                          sep = "")
directory.params <- paste(directory.models,
                          'NFkB/Input/parameters_names.csv',
                          sep = "")

dir.create(path = directory.output, recursive = TRUE, showWarnings = FALSE)

SM <- as.matrix(read.csv(directory.SM, sep = "\t", header = FALSE))
params <- read.csv(directory.params, header = FALSE, sep = "\t") ## HORIZONTAL

if( ncol(SM) != nrow(params)){
  print("Error. Wrong number of parameters")
}

x <- SMC(S = SM, 
         labels = t(params),
         names = t(params),
         zeta = zeta,
         delta = delta,
         dir.output =directory.output )


#### JAKSTATGB EXAMPLE ####
directory.folder <- paste(directory.models,
                          'JAKSTATGB/',
                          sep = "")
directory.SM     <- paste(directory.folder,
                          'Output/Experiment_Jaccobo/sense.csv',
                          sep = "")

directory.output <- paste(directory.folder,
                          'Output/Experiment_Jaccobo/',
                          "clustering-2/",
                          'delta/',
                          '/', 
                          sep = "")
directory.params <- paste(directory.folder,
                          'Input/parameters_names.csv',
                          sep = "")

dir.create(path = directory.output, recursive = TRUE, showWarnings = FALSE)

SM <- as.matrix(read.csv(directory.SM, sep = "\t", header = FALSE))
params <- read.csv(directory.params, header = FALSE, sep = "\t") ## HORIZONTAL

if( ncol(SM) != nrow(params)){
  print("Error. Wrong number of parameters")
}

x <- SMC(S = SM, 
  labels = t(params),
  names = t(params),
  zeta = zeta,
  delta = delta,
  dir.output =directory.output )


#### NFKB EXAMPLE ####
directory.folder <- paste(directory.models,
                          'tjetka/Output/',
                          sep = "")
directory.SM     <- paste(directory.folder, 'sense.csv', sep = "")
directory.FIM     <- paste(directory.folder, 'FIM.csv', sep = "")
directory.output.SM <- paste(directory.folder,
                          "clustering/",
                          'SM',
                          '/', 
                          sep = "")
directory.output.FIM <- paste(directory.folder,
                          "clustering/",
                          'FIM',
                          '/', 
                          sep = "")

### case SM
dir.create(path = directory.output.SM, recursive = TRUE, showWarnings =   FALSE)
SM <- as.matrix(read.csv(directory.SM, sep = "\t", header = FALSE))

x.SM <- SMC(S = SM, 
  zeta = zeta,
  delta = delta,
  dir.output =directory.output.SM )

directory.output <- directory.output.SM
x <- x.SM

### case FIM
dir.create(path = directory.output.FIM, recursive = TRUE, showWarnings = FALSE)
FIM <- as.matrix(read.csv(directory.FIM, sep = "\t", header = FALSE))

x.FIM <- SMC(FIM = FIM, 
            zeta = zeta,
            delta = delta,
            dir.output =directory.output.FIM )

directory.output <- directory.output.FIM
x <- x.FIM


#### MAIN CODE ####

cluster <- clusterident.SMC(x)

directory.dendrogram <- paste(directory.output, "_dendrogram", ".pdf", sep = "")
pdf(directory.dendrogram, width = 12, height = 6)
plotDendogram.SMC(x,
                  cluster = cluster,
                  fig      = c(0,1,0,1),
                  new      = FALSE,
                  labels.args = list(cex = 0.5)
)
dev.off()

directory.fi <- paste(directory.output, "_fi", ".pdf", sep = "")
pdf(directory.fi, width = 6, height = 6)
p <- plotFI.SMC(x,
                cluster = cluster,
                fig      = c(0,1,0,1),
                new      = FALSE
)
dev.off()

directory.barplot <- paste(directory.output, "_barplot", ".pdf", sep = "")
pdf(directory.barplot, width = 12, height = 6)
barplot.SMC(x, cluster = cluster, fig = c(0,1,0,1))
dev.off()

directory.barplot <- paste(directory.output, "_barplot", ".pdf", sep = "")
pdf(directory.barplot, width = 12, height = 6)
barplot.SMC(x, cluster = cluster, fig = c(0,1,0,1))
dev.off()

directory.cov <- paste(directory.output, "_cov", ".pdf", sep = "")
pdf(directory.cov, width = 6, height = 6)
plotEV.Cov.SMC(x, cluster = cluster, fig = c(0,1,0,1))
dev.off()

directory.cor <- paste(directory.output, "_cor", ".pdf", sep = "")
pdf(directory.cor, width = 6, height = 6)
plotEV.Cor.SMC(x, cluster = cluster, fig = c(0,1,0,1))
dev.off()

directory.sa <- paste(directory.output, "_sa", ".csv", sep = "")
data.sa <- sa.SMC(x, cluster = cluster, header=FALSE)
write.csv(x = data.sa, file = directory.sa)
