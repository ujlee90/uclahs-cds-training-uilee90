#To this point all the analyses have been performed using standard R functions.  
#Fortunately the BoutrosLab repository has a variety of functions that will aid in making the plots and statistics.  
#Modify your code to use:
library(BoutrosLab.statistics.general)
library(BoutrosLab.plotting.general)
  #a) BoutrosLab.statistics.general	for all p-value extraction functions

read.txt <- function(x, sep){
  input <- read.table(x, 
                      header=TRUE, 
                      sep=sep,
                      na.strings = "NA",
                      stringsAsFactors = default.stringsAsFactors());
}
input1 <- read.txt("input1.txt","\t")
input2 <- read.txt("input2.txt","\t")
head(input1)
head(input2)

# sort the files
data.sort <- function(x, split, sortID=GeneID){
  tmp.data <- t(matrix(
    unlist(strsplit(as.character(x$GeneID),
                    split=split))));
  geneid <- order(as.numeric(tmp.data));
  x.sort <- x[geneid,];
  return(x.sort);
}
input1.sort <- data.sort(x=input1, split="_at")
input2.sort <- data.sort(x=input2, split="_at")
input.all.merge <- merge(input1.sort, input2.sort, key="GeneID")
input.all.sort <- data.sort(x=input.all.merge, split="_at")
head(input.all.sort,25)


set.seed(12345)

# t.test p.values with BoutrosLab.statistics.general
t.test.p.values <- apply(input.all.sort[,2:13], 1, function(x)
  BoutrosLab.statistics.general::get.ttest.p(
    as.numeric(x),
    group1 = c(1:3),
    group2 = c(4:12)
  )
)
head(t.test.p.values)

# Wilcox test p.values with BoutrosLab.statistics.general
wilcox.test.p.values <- apply(input.all.sort[,2:13], 1, function(x)
  BoutrosLab.statistics.general::get.utest.p(
    as.numeric(x), 
    group1 = c(1:3),
    group2 = c(4:12)
  )
)
head(wilcox.test.p.values)

# fold change with BoutrosLab.statistics.general
foldchange.values <- apply(input.all.sort[,2:13], 1, function(x)
  BoutrosLab.statistics.general::get.foldchange(
    as.numeric(x), 
    group1 = c(1:3),
    group2 = c(4:12),
    logged = FALSE
  )
)
head(foldchange.values)

#b) BoutrosLab.plotting.general	for all plots

# Histogram of t.test p.values with BoutrosLab.plotting.general
BoutrosLab.plotting.general::create.histogram(
  x = as.numeric(t.test.p.values),
  main = "Student's t test p-values of 500 genes for tumor type A vs B",
  main.cex = 1.2,
  xlab.label = "p-values",
  xlab.cex = 1.25,
  xlimits = c(0,1),
  ylimits = c(0, 15),
  ylab.cex = 1.25,
  ylab.label = "Frequency",
  type = "count",
  xat =  c(seq(0, 1, 0.2)),
  breaks = seq(0, 1, 0.01),
)

# Histogram of wilcox p.values with BoutrosLab.plotting.general
BoutrosLab.plotting.general::create.histogram(
  x = as.numeric(wilcox.test.p.values),
  main = "Wilcox t test p-values of 500 genes for tumor type A vs B",
  main.cex = 1.2,
  xlab.label = "p-values",
  xlab.cex = 1.25,
  xlimits = c(0,1),
  ylimits = c(0, 80),
  ylab.cex = 1.25,
  ylab.label = "Frequency",
  type = "count",
  xat =  c(seq(0, 1, 0.2)),
  breaks = seq(0, 1, 0.025),
)
  

# Histogram of t.test p.values with BoutrosLab.plotting.general
BoutrosLab.plotting.general::create.histogram(
  x = as.numeric(foldchange.values),
  main = "Fold change of 500 genes for tumor type A vs B",
  main.cex = 1.2,
  xlab.label = "Fold change",
  xlab.cex = 1.25,
  ylab.cex = 1.25,
  ylab.label = "Frequency",
  type = "count",
  breaks = 100,
)

#Make these modifications for all questions from Q2 onwards.