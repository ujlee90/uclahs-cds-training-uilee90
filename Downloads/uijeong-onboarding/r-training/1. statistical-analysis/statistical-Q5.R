#As noted above it is often hard to identify an appropriate analytical statistical test. 
#This may be true because assumptions are violated, or because no standard test exists. 
#In either case, a common solution is to use permutation-based approaches, like a bootstrap.  
#Your final task is to develop a bootstrap-like test in R.  This may be quite challenging at first.
library(BoutrosLab.plotting.general)
  
read.txt <- function(x, sep){
  input <- read.table(x, 
                       header=TRUE, 
                       sep=sep,
                       na.strings = "NA",
                       stringsAsFactors = default.stringsAsFactors())
}
input1 <- read.txt("input1.txt","\t")
input2 <- read.txt("input2.txt","\t")

input1_new <- input1
input2_new <- input2
input_merge <- merge(input1_new, input2_new, key="GeneID")
head(input_merge)

#1. Calculate the median of the first three columns for each gene

median.first.3 <- apply(input_merge, 1, function(x) median(as.numeric(x[2:4])))

head(median.first.3)

#2. Use a permutation test to estimate the expected value for each gene:

 
bootstrap.larger <- NULL
bootstrap.smaller <- NULL
bootstrap <- NULL
num = 1
while (num < 1001) {           #d. Repeat a.-c. 1000 times
  median.first.3 <- apply(input_merge, 1, function(x) median(as.numeric(x[2:4])));
  random.3 <- apply(input_merge, 1, function(x) sample(as.numeric(x[2:13]),3,replace = TRUE)); #a. Randomly select three columns from amongst all 12
  median.random <- apply(random.3, 2, function(x) median(x)); #b. Calculate their median
  median.first3.random <- cbind(median.first.3, median.random)
  bootstrap.larger <- apply(median.first3.random, 1, function(x) x[2] > x[1]) #c. Determine if this value is larger or smaller than that of the first 3 columns
  bootstrap.smaller <- apply(median.first3.random, 1, function(x) x[2] < x[1])
  bootstrap <- cbind(bootstrap, bootstrap.larger);
  num <- num + 1
}
head(median.first3.random)
head(bootstrap)
frequency <- NULL
frequency <- apply(bootstrap, 1, function(x) sum(x))
frequency

#3. Use the frequencies in 2. to estimate a p-value for each gene 

estimate.frequency <- frequency/1000
estimate.frequency

#4. Perform a false-discovery adjustment on the p-values (?p.adjust)

estimate.frequency.fdr <- p.adjust(estimate.frequency, method = "fdr", n=1000)

#5. Write your results (gene ID, observed median, expected median, p-value, adjusted p-value) to file in a tab-delimited format

total.results <- data.frame(GeneID=input_merge$GeneID, observed.median= median.random, expected.median=median.first.3, p.value=estimate.frequency, adjusted.p.value=estimate.frequency.fdr)
head(total.results)

#6. Plot a histogram of the (unadjusted) p-values. What does this tell you?
hist(total.results$p.value,
     main = "Bootstrap p-value",
     freq = TRUE,
     xlab = "p-values",
     breaks = seq(0,1, 0.01))

########Q6.BoutrosLab.plotting.general	for all plots
# Histogram of Bootstrap p.values with BoutrosLab.plotting.general
BoutrosLab.plotting.general::create.histogram(
  x = as.numeric(total.resutls$p.value),
  main = "Bootstrap p-values",
  main.cex = 1.2,
  xlab.label = "p-value",
  xlab.cex = 1.25,
  ylab.label = "Frequency",
  ylab.cex = 1.25,
  type = "count",
  xat =  c(seq(0, 1, 0.2)),
  breaks = 100,
)

#What does this tell you?
  # > The Bootstrap randomly selects samples with replacement from a small population, and it forms a normal distribution and obtains a confidence interval(95%). 
  #   Above the results, there is the probability that the median of the randomly three selected values (from the 12 types) is greater than the median of the the three samples (tumor type A), which is obtained by repeating it 1000 times for each gene.
