#Okay, your next question might be if a t-test was inappropriate for this analysis.  
#Repeat the above comparison using a Wilcoxon test and fold-changes.  
#Create suitable plots to compare the results and write a paragraph describing the similarities/differences.

# open the files
read.txt <- function(x, sep){
  input1 <- read.table(x, 
                       header=TRUE, 
                       sep=sep,
                       na.strings = "NA",
                       stringsAsFactors = default.stringsAsFactors())
}
read.txt("input1.txt","\t")
read.txt("input2.txt","\t")

# sort the files
data.sort <- function(x, split, sortID=GeneID){
  tmp.data <- t(matrix(
    unlist(strsplit(as.character(x$GeneID),
    split=split))));
  geneid <- order(as.numeric(tmp.data));
  x.sort <- x[geneid,];
  return(x.sort);
}
input1.sort <- data.sort(x=input1.new, split="_at")
head(input1.sort)
input2.sort <- data.sort(x=input2.new, split="_at")
input.all.merge <- merge(input1.sort, input2.sort, key="GeneID")
input.all.sort <- data.sort(x=input.all.merge, split="_at")
head(input.all.sort,25)
sum(is.na(input.all.sort))

#wilcox.test and fold-changes

wilcox.fold <- NULL

for ( i in 1:nrow(input.all.sort)) {
  x <- wilcox.test(
    as.numeric(input.all.sort[i,c(2:4)]), 
    as.numeric(input.all.sort[i,c(5:13)]),
    alternative="two.sided",
    paired = FALSE,
    conf.int = FALSE,
    conf.level = 0.95,
    exact = FALSE
    )$p.value;
  genename <- as.vector(input.all.sort[i,1]);
  fold.changes <- as.numeric(mean(unlist(input.all.sort[i,c(2:4)]))) / as.numeric(mean(unlist(input.all.sort[i,c(5:13)])))
  fold.changes <- log2(unlist(fold.changes))
  wilcox.fold <- rbind(wilcox.fold, c(genename, x, fold.changes)
  );
}

head(wilcox.fold)
wilcox.fold <- data.frame(Gene_ID = wilcox.fold[,1], wilcox.p.value = as.numeric(wilcox.fold[,2]), fold.change = as.numeric(wilcox.fold[,3])) #change to dataframe
head(wilcox.fold)

#histogram of the wilcox.p.value
hist(wilcox.fold$wilcox.p.value,
     main = "Wilcoxon test p-values of 500 genes: Tumor type A vs B",
     freq = TRUE,
     xlab = "p-values",
     breaks = seq(0,1, 0.1))

#histogram of the fold change
hist(wilcox.fold$fold.change,
     main = "Log2 fold change of 500 genes: Tumor type A vs B",
     freq = TRUE,
     xlab = "Log2 fold change",
     breaks = 100)


library(ggplot2)
volcano = ggplot(data = wilcox.fold, aes(x = fold.change, y = -1*log10(wilcox.p.value)))
volcano + geom_point()

#>There are the results between group1[first three column] and group2[last nine column] obtained by t-test and Wilcoxon test for each gene. 
#> However, I think that 'n' is not enough to know this data meaningful through the value of p-value for each gene. With some graphs of the Wilcoxon test, there are many p-values over 0.05 and the results say that the score of the two groups is not different for the gene.
#> In addition, it is known that the fold change/log2-fold change shows gene expression and is used with DEG analysis. The number shown in the input files does not know what the accuracy means, but . It is visually seen as a scatter plot or a volcano plot.



#To this point you've been doing these analyses using a loop (for or while statements).  
#You'll next want to learn to use the R function "apply" to generate a vector of p-values (t- and u-tests) and of fold-changes.  
#Remember, the vectors should only contain p-values, not any other information.  
#You will need to figure out how to do this using apply, and this will require you to create what's called a "wrapper function".


ttest.p.value <- apply(input.all.sort[,2:13], 1, function(x){
  tmp <- as.vector(x);
  p.value <-t.test(
    as.vector(tmp[c(1:3)]), 
    as.vector(tmp[c(4:12)]),
    alternative="two.sided",
    paired = FALSE,
    conf.int = TRUE,
    conf.level = 0.95,
    exact = TRUE
    )$p.value;
  wilcox <- wilcox.test(
    as.numeric(tmp[c(1:3)]), 
    as.numeric(tmp[c(4:12)]), 
    alternative="two.sided",
    paired = FALSE,
    conf.int = FALSE,
    conf.level = 0.95,
    exact = FALSE
    )$p.value;
  fold.changes <- mean(as.numeric(tmp[c(4:12)]), na.rm=TRUE) / mean(as.numeric(tmp[c(1:3)]), na.rm=TRUE);
  return(c(p.value, wilcox));
  }
)
ttest.p.value <- t(ttest.p.value)
head(ttest.p.value)
