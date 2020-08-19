#There are a number of statistical challenges involved in executing so many statistical tests.  
#The primary one is called the "multiple testing" problem.  What is it? 
  # > In a hypothesis test, there is a bogus significant result(about 5%). If thousands of tests are run, the number of false alarms increases dramatically. 
  #   For instance, If I use the standard alpha level of 5% in 10,000 separate hypothesis tests, I would get around 500 significant results, most of which will be false alarms. 
  #   When I run multiple hypothesis tests, this large number of false alarms produced, and it is called the multiple testing problem. (Or multiple comparisons problem). significant result. 

#Find R commands to adjust for multiple-testing using the two best-known adjustments (FDR and Bonferroni).
  # > p.adjust(data$p.value, 'fdr'), p.adjust(data$p.value, 'bonferroni')
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

p.values <- apply(input.all.sort[,2:13], 1, function(x){
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
p.values <- t(p.values)
head(p.values)
  
fdr.ttest <- p.adjust(p.values[,1], method = 'fdr', n=500)
fdr.wilcox <- p.adjust(p.values[,2], method = 'fdr', n=500)

bonferroni.ttest <- p.adjust(p.values[,1], method = 'bonferroni', n=500)
bonferroni.wilcox <- p.adjust(p.values[,2], method = 'bonferroni', n=500)

multiple.testing <- data.frame(ttest.p.value= as.numeric(p.values[,1]), 
                            wilcox.p.vlaue= as.numeric(p.values[,2]),
                            fdr.ttest=as.numeric(fdr.ttest), 
                            fdr.wilcox=as.numeric(fdr.wilcox),
                            bonferroni.ttest=as.numeric(bonferroni.ttest),
                            bonferroni.wilcox=as.numeric(bonferroni.wilcox))
head(multiple.testing)
table(multiple.testing)

#Recreate the histograms above.  What does this tell you?
hist(multiple.testing$fdr.ttest,
     main = "Unpaired t-test FDR adjusted p-value of 500 genes",
     xlab = "p-values",
     breaks = seq(0,1,0.05)
     )
hist(multiple.testing$fdr.wilcox,
     main = "Wilcoxon test FDR adjusted p-value of 500 genes",
     xlab = "p-values",
     breaks = seq(0,1,0.05)
     )
hist(multiple.testing$bonferroni.ttest,
     main = "Unpaired t-test Bonferroni adjusted p-value of 500 genes",
     xlab = "p-values",
     breaks = seq(0,1,0.05)
     )
hist(multiple.testing$bonferroni.wilcox,
     main = "Wilcoxon test Bonferroni adjusted p-value of 500 genes",
     xlab = "p-values",
     breaks = seq(0,1,0.05))

#With the FDR and Bonferroni corrections, It is corrected for false positives that arise from multiple testing with inference statistics. 
#Above the results, I can see that the FDR correction is more tolerant than the Bonferroni about a t-test and a Wilcoxon test.
#A number of genes with p-values(t-test) between 0.1 - 0.15 are located when performing the FDR correction , and it could be statistically significant. 
#On the other hand, the Bonferroni correction is more stringent and corrects most p-values for the t-test and Wilcoxon test to 1. 
#Thus, there would be no genes which have have any statistical significance by using the Bonferroni correction.

  
#Hopefully you have been saving your code in scripts as suggested in the bl-manulal. 
#The code you create is an important resource to you and your colleagues and it is for this reason that we want to be able to track changes to this code and share it with others. 
#One way to do this is via GitHub. Please read and follow the steps in the git/GitHub guide on Box https://uclahs.box.com/s/qv3olx8qnkuz4znl1rt6dkah4erbqce8
#You will want to add your code and output from questions 1-4 to your local r-statistics folder within your <username>-training repository before pushing to Github.