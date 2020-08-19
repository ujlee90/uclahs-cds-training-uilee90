install.packages("BoutrosLab.utilities_1.9.10.tar.gz", type="source", dependencies=TRUE, repos=NULL)
install.packages("BoutrosLab.statistics.general_2.1.3.tar.gz", type="source", dependencies=TRUE, repos=NULL)

#Read the file AHR-test-file.txt
AHRtest <- read.table(file = "AHR-test-file.txt",
                    header = TRUE, 
                    sep = "\t",
                    na.strings = "NA",
                    stringsAsFactors = default.stringsAsFactors())
AHRtest

AHR.test <- AHRtest # copy the content in the file
sum(is.na(AHR.test)) # check the "NA" 

#perform a t-test between control and treated

var.test(AHR.test$Control, AHR.test$Treated, na.rm=TRUE)$p.value # p-value > 0.05
t.test(AHR.test$Control, AHR.test$Treated, alternative = ="two.sided", conf=0.95, var.equal=FALSE, paired=FALSE)

#perform a wilcoxon test between control and treated

wilcox.test(AHR.test$Control, AHR.test$Treated, alternative="two.sided", paired=FALSE, conf.level=0.95, exact=FALSE)

#calculate a fold-change between control and treated

fold.change <- log2(mean(na.omit(AHR.test$Treated))/mean(na.omit(AHR.test$Control)))
fold.change


###Theory
#Different type of t-tests
  # >indepentent two sample t-test
  # >paired(dependent) sample t-test
  # >one sample t-test
  
#Different two-sample type of tests
  # >Parametric Test(normal distribution) -> t-test(independent, paired)
  # >Non-parametric Test(non-normal distribution or unknown distribution) -> wilcoxon rank sum test & mann-whitney U-test, wilcoxon signed rank test

#Which test to use when?
  # >paired t-test : use this when the same subjects participate in both conditions of the experiment.
  # >independent t-test : use this when you have two different groups of subjects, one group performing one condition in the experiment, and the other group performing the other condition.
  # >wilcoxon test : use this when you have two groups with non-normal distribution or unknown distribution
  # >ANOVA test : use this when you have over three different groups of subjects  

###Coding
  #Good file-opening semantics : The biggest important thing is that files and their contents are fully opened without missing ones or errors.
  #                              We have to confirm the file pathway using getwd(),setwd(), and the file format(read.csv, read.table, readxl::read_excel, etc).
  #                              The identification of the contents in the file is also needed because of some options such as header=TRUE, sep="", col.names=c(), na.string="", etc.
  #R help : >help()
  #Prameterization in R functions : Any parameters are set as arguments to the function. The code wraps all of specified arguments into a function.
  #                                 This is exemplified above; function(x, ...). 








