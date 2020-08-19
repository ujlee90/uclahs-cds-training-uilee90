#In bioinformatics it is often hard to identify an appropriate analytical statistical test. This may be true because assumptions are violated, or because no standard test exists. In either case, a common solution is to use permutation-based approaches, like a bootstrap.

#In this problem you will be comparing a number of different statistical approaches for comparing two groups.  The biological question here involves two tumour subtypes, A and B.  You have three samples of subtype A and nine samples of subtype B.  You have used a small microarray platform to measure the mRNA levels of 500 genes on each of these 12 tumours.  You wonder which genes differ between the two tumour subtypes.

#Your first input data file (input1.txt) contains one column for the gene identifier and three columns of numeric data. 
#Your second input data file (input2.txt) contains one column for the gene identifier and nine columns of numeric data. 
#For both files, each column represents a biological sample (human tumour) and each row represents a gene.  
#The three numeric columns in input1.txt represent tumour-type A and the nine numeric columns in input2.txt represent tumour-type B.  
#In statistical terms, you need to determine if the three columns in input1.txt represent a random sample from the overall data for each gene.

#Your first steps to answer this question are:
#1. Read the two input files
read.txt <- function(x, sep){
  input <- read.table(x, 
                       header=TRUE, 
                       sep=sep,
                       na.strings = "NA",
                       stringsAsFactors = default.stringsAsFactors());
}
input1 <- read.txt("input1.txt","\t");
input2 <- read.txt("input2.txt","\t");
head(input1)
head(input2)

#2. Combine the two files into one file that contains the data for all 12 tumours. Make sure that the three columns in input1.txt precede the nine columns in input2.txt. Do this in the following two ways, and verify that they produce the same result:
#a) Sort each file individually, and then use the cbind function
install.packages("dplyr");
library(dplyr);

input1.new <- input1
input2.new <- input2

data.sort <- function(x, split, sortID=GeneID){
          tmp.data <- t(matrix(
            unlist(strsplit(as.character(x$GeneID),
              split=split))));
          geneid <- order(as.numeric(tmp.data));
          x.sort <- x[geneid,];
          return(x.sort);
}
input1.sort <- data.sort(x=input1.new, split="_at");
input2.sort <- data.sort(x=input2.new, split="_at");

input.total <- cbind(input1.sort[,1:4], input2.sort[2:10]);
head(input.total,25)

#b) Use only the merge function
input.all.merge <- merge(input1.sort, input2.sort, key="GeneID");
input.all.sort <- data.sort(x=input.all.merge, split="_at");
head(input.all.sort,25)
sum(is.na(input.all.sort))


#3. Perform a t-test comparing the first three tumours to the last nine tumours for *each* gene using a for-loop

gene.p.value <- NULL;

for ( i in 1:nrow(input.all.sort)) {
    x <- t.test(
      as.vector(input.all.sort[i,c(2:4)]), 
      as.vector(input.all.sort[i,c(5:13)]),
      alternative="two.sided",
      paired = FALSE,
      conf.int = TRUE,
      conf.level = 0.95,
      exact = TRUE
      )$p.value;
    genename <- as.vector(input.all.sort[i,1]);
    gene.p.value <- rbind(gene.p.value, c(genename, x)
  );
}
#colnames(gene.p.value) <- c('Gene_ID', 'P.value')

gene.p.value <- data.frame(
  Gene_ID = gene.p.value[,1], 
  P.value = as.numeric(gene.p.value[,2])); #change to dataframe
head(gene.p.value)

shapiro.test(gene.p.value$P.value) #confirm the normal distribution of p.value

#4. Plot a histogram of the p-values
#method 1.
str(gene.p.value)
hist(gene.p.value$P.value,
     main = "Student's t test p-value of 500 genes",
     freq = TRUE,
     xlab = "p-values",
     breaks = seq(0,1, 0.01));

#method 2.
#install.packages("ggplot2")
library(ggplot2)
ggplot(data=gene.p.value, 
       aes(x=P.value)) + geom_histogram();
ggsave("gene.p.value.png")


#5. Are your axis labels rotated 90 degrees?  If so, fix this.
#method 1.
hist(gene.p.value$P.value,
     main = "Student's t test p-value of 500 genes",
     freq = TRUE,
     xlab = "p-values",
     breaks = seq(0,1, 0.01),
     las = 2
     );

#method 2.
ggplot(data = gene.p.value, aes(x = P.value)) + geom_histogram() +
  theme(axis.text.x=element_text(angle = 90));


#6. Your histogram might look a bit weird in normal space, consider plotting it in log-space

log <- -log(as.numeric(gene.p.value[,2]),base=10);
log.gene <- cbind('Gene_ID'=gene.p.value[,1],'log.p.value'=log);
log.p.value <- data.frame(Gene_ID = log.gene[,1], Log.p.value = as.numeric(log.gene[,2]));
hist(log.p.value$Log.p.value);

#7. What does this distribution tell you?
  # > The p-value of shapiro-wilk normality test is lower than the significance level alpha = 0.05. It means that this is not normally distributed.  
  #   In the histogram, the distribution of frequency of p.value is fairly evenly distributed from 0 to 1. 
  #   This rejects the null hypothesis(a variable is normally distributed in some population) and follow the alternative hypothesis.
  #   It could be concluded that it is needed to conduct another test such as wilcoxon test.

  