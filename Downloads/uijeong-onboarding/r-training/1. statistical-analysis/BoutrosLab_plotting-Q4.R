#WORKING WITH LARGE DATASETS: BoutrosLab.plotting.general
#=======
#  4.
#- Take a look at "Q4_SampleOutput.tiff"
#- Re-create this figure using BoutrosLab.plotting.general functions. 
#- The data file you will use is named "Q4_HetStudy_data.txt"
library(BoutrosLab.plotting.general);
library(dplyr);

HetStudy_data <- read.delim("Q4/Q4_HetStudy_data.txt", header = TRUE, sep ="\t");
head(HetStudy_data)

#Understanding the figure: (more information on the project can be found by reading the paper https://www.nature.com/articles/ng.3315)
#- The fraction plot at the top of the figure indicates the proportion of mutations per bin. (Note: your results will not look identical to the plot provided).
#- The literature covariate bars labelled "Baca", "Berger", and "Weischenfeldt" are studies. The colours indicate whether or not a mutation in that window of the genome was seen in that study.
#- The large plot of small coloured bars indicates locations in chromosomes (x-axis) in different samples (y-axis) of mutations. CTX = chromosomal translocation, ITX = intrachrosomosomal translocation, INV = inversion.
#- The right-side covariates identify clinical qualities of the samples. The first covariate indicates which patient the samples come from. The second covariate indicates whether the sample was collected by surgery or biopsy. The third covariate indicates the Gleason score of the sample. The fourth covariate indicates how the tissue samples were preserved. 

#Understanding the data:
#  - The patient IDs are identified in the column headers.
#- The cohorts are such that samples "CPCG0001", "CPCG0003", "CPCG0006", "CPCG0020", and "CPCG0042" are Bx (biopsy), and all others are Sx (surgery)
#- For tissue type, all "F0" samples are Frozen, while all other samples are FPPE
#- The Gleason score can be matched to the sample as shown in the sample output plot
#- The literature covariate data (Baca, Berger, and Weischenfeldt) are listed in the last three columns of the data file.
#- The fraction plot at the top is found by the proportion of mutations per patient
#- The notation for None, CTX, ITX and INV is determined by the numeric coding of 0, 1, 2, 3.

#Colour schemes used:
#  -main plot: c('white', 'cornflowerblue','darkolivegreen4','darkred')
#-Cohort: c('royalblue', 'pink')
#-Gleason score: c('3+4' = 'yellow', '4+3' = 'orange','4+4' = 'red')
#-Tissue type: c('Frozen' = colours()[532],'FFPE' = colours()[557])
#-Patient ID: c( 
  'CPCG0001' = 'blue',
  'CPCG0003' = 'purple',
  'CPCG0006' = 'green',
  'CPCG0020' = 'orange',
  'CPCG0042' = 'yellow',
  'CPCG0099' = 'black',
  'CPCG0102' = 'wheat4',
  'CPCG0103' = 'green4',
  'CPCG0183' = 'grey',
  'CPCG0184' = 'red4')
#-Literature: you may choose your own appropriate colour scheme

#Sample Order for plotting a heap map:
Q4_main_data.mtx <-select(Q4_HetStudy_data, CPCG0001F0:CPCG0184P4);
Q4_BBW_data <-select(Q4_HetStudy_data, Baca:Weischenfeldt);
head(Q4_BBW_data)
sample.order <- c(
    "CPCG0001F0",
    "CPCG0003F0",
    "CPCG0006F0",
    "CPCG0020F0",
    "CPCG0042F0",
    "CPCG0099F0",
    "CPCG0099P1",
    "CPCG0102F0",
    "CPCG0102P1",
    "CPCG0102P2",
    "CPCG0103P7",
    "CPCG0103F0",
    "CPCG0103P2",
    "CPCG0103P1",
    "CPCG0103P4",
    "CPCG0103P3",
    "CPCG0103P8",
    "CPCG0103P5",
    "CPCG0103P6",
    "CPCG0183F0",
    "CPCG0183P2",
    "CPCG0183P1",
    "CPCG0183P3",
    "CPCG0184P3",
    "CPCG0184P1",
    "CPCG0184F0",
    "CPCG0184P2",
    "CPCG0184P4"
  );

#############################################################
#main heatmap
sampleorder.data.mtx <- Q4_main_data.mtx[,sample.order];
head(sampleorder.data.mtx)
str(sampleorder.data.mtx)

col.line.v <- NULL
i = 1
while (i < 23) {
  name.temp <- c("chr");
  name.temp2 <- paste(name.temp, as.character(i),":", sep = "");
  temp <-length(grep(name.temp2, row.names(sampleorder.data.mtx)));
  col.line.v <- c(col.line.v, temp);
  i = i + 1              
}

col.line.v <- c(col.line.v, length(grep("chrX", row.names(sampleorder.data.mtx))), length(grep("chrY", row.names(sampleorder.data.mtx))));
col.line.v
cum.col.line.v <- c(cumsum(col.line.v));
cum.col.line.v
xlabel.v.1 <- cum.col.line.v - col.line.v/2;
str(xlabel.v.1)

#making the heatmap
mainQ4_heatmap <- create.heatmap(
  x = sampleorder.data.mtx,
  clustering.method = "none",
  scale.data = FALSE,
  colour.scheme = c('white', 'cornflowerblue','darkolivegreen4','darkred'),
  total.colours = 5,
  xaxis.lab = as.vector(c(c(1:22, "X", "Y"))),
  xaxis.cex = 0.75,
  yaxis.lab = (c(substring(sample.order, 9, 10))),
  yaxis.cex = 0.75,
  xat = xlabel.v.1,
  xaxis.rot = 0,
  print.colour.key = FALSE,
  covariates.col.lines = TRUE,
  grid.col = TRUE,
  grid.row = TRUE,
  col.lwd =1.25,
  row.lwd = 1.25,
  col.colour = 'black',
  col.lines = cumsum(col.line.v),
  row.lines = c(1, 2, 3, 4, 5, 7, 10, 19, 23) +0.5,
  force.grid.col = TRUE,
)

mainQ4_heatmap

#############################################################
#creading the bar plot on the top

head(Q4_main_data.mtx)
bar_plot_data.v <- apply(Q4_main_data.mtx, 1, function(x) (28-sum((x == 0 )))/ 28 );
bar_plot_data.v

bar_plot_data.df <- data.frame( 
  y = bar_plot_data.v,
  x = c(1:length(bar_plot_data.v))
);
bar_plot_data.df

#making the barplot
Q4_bar_plot <- create.barplot(
  formula = y ~ x,
  data = bar_plot_data.df,
  ylab.label =  "Fraction",
  xlab.label = NULL,
  xaxis.cex = 0,
  yaxis.cex = 2,
  yat = seq(0, 1.0, 0.5),
  ylimits = c(0, 1),
  xlimits = c(0:length(bar_plot_data.v)),
  xaxis.tck = 0,
  xat = NULL,
  xaxis.lab = rep("",length(bar_plot_data.v))
)

Q4_bar_plot

#############################################################
#creating covariate_1(Baca, Berger, and Weischenfeldt)

covariate_1.data <-cbind(Q4_BBW_data[,1], Q4_BBW_data[,3], Q4_BBW_data[,2]);
head(covariate_1.data)

#colors
covaraite_1.colours <- c('white','red','deepskyblue3','mediumseagreen');

covariate_1.data[covariate_1.data == 2] <- 4;
covariate_1.data[covariate_1.data == 1] <- 2;
covariate_1.data[covariate_1.data == 0] <- 1;
head(covariate_1.data)

#making the heatmap
covariate_1 <- create.heatmap(
  x = covariate_1.data,
  clustering.method = 'none',
  colour.scheme = as.vector(covaraite_1.colours),
  total.colours = 5,
  grid.row = TRUE,
  grid.col = TRUE,
  col.lwd =1.25,
  row.lwd = 1.25,
  col.colour = 'black',
  row.colour = 'black',
  yaxis.tck = 0,
  print.colour.key = FALSE,
  force.grid.col = TRUE,
  col.lines = cum.col.line_vector
);

covariate_1

#####################################################
#creating covariate_2(None, CTX, ITX, INY)

covariate_2 <- create.heatmap(
  x = data.frame(1,2,3,4),
  clustering.method = 'none',
  colour.scheme = c('white', 'cornflowerblue','darkolivegreen4','darkred'),
  total.colours = 5,
  xaxis.lab = c('None', 'CTX', 'ITX', 'INV'),
  grid.col = TRUE,
  col.lwd =1.25,
  col.colour = 'black',
  xaxis.rot = 0,
  xaxis.cex = 1.5,
  yaxis.cex = 0,
  same.as.matrix = TRUE,
  print.colour.key = FALSE,
)

covariate_2


#############################################################
#creating covariate_3 on the right side

HetStudy_clinical.data <- read.delim("Q4/Q4_HetStudy_clinical.txt", header = TRUE);
head(HetStudy_clinical.data)

covariate_3.data <- as.matrix(
  cbind(
    substring(sample.order, 1, 8),
    HetStudy_clinical.data[,c(2:4)]
  )
);

head(covariate_3.data)
str(covariate_3.data)

covaraite_3.colours <- c('royalblue', 'pink', 'yellow', 'orange', 'red', 'blue', 'purple', 'green', 'black', 'wheat4', 'green4', 'grey', 'red4', 'bisque2', 'bisque4');

#the heatmap expects numeric data
covariate_3.data[covariate_3.data == 'CPCG0001'] <- 6;
covariate_3.data[covariate_3.data == 'CPCG0003'] <- 7;
covariate_3.data[covariate_3.data == 'CPCG0006'] <- 8;
covariate_3.data[covariate_3.data == 'CPCG0020'] <- 4;
covariate_3.data[covariate_3.data == 'CPCG0042'] <- 3;
covariate_3.data[covariate_3.data == 'CPCG0099'] <- 9;
covariate_3.data[covariate_3.data == 'CPCG0102'] <- 10;
covariate_3.data[covariate_3.data == 'CPCG0103'] <- 11;
covariate_3.data[covariate_3.data == 'CPCG0183'] <- 12;
covariate_3.data[covariate_3.data == 'CPCG0184'] <- 13;
covariate_3.data[covariate_3.data == 'Sx'] <- 2;
covariate_3.data[covariate_3.data == 'Bx'] <- 1;
covariate_3.data[covariate_3.data == '3+4'] <- 3;
covariate_3.data[covariate_3.data == '4+3'] <- 4;
covariate_3.data[covariate_3.data == '4+4'] <- 5;
covariate_3.data[covariate_3.data == 'Frozen'] <- 14;
covariate_3.data[covariate_3.data == 'FFPE'] <- 15;
covariate_3.data <-as.data.frame(covariate_3.data);

covariate_3.data[,1] <- as.numeric(as.character(covariate_3.data[,1]));
covariate_3.data[,2] <- as.numeric(as.character(covariate_3.data[,2]));
covariate_3.data[,3] <- as.numeric(as.character(covariate_3.data[,3]));
covariate_3.data[,4] <- as.numeric(as.character(covariate_3.data[,4]));

head(covariate_3.data)
covariate_3 <- create.heatmap(
  x = t(covariate_3.data),
  same.as.matrix = FALSE,
  clustering.method = 'none',
  colour.scheme = as.vector(covaraite_3.colours),
  yaxis.tck = 0,
  total.colours = 16,
  grid.row = TRUE,
  grid.col = TRUE,
  print.colour.key = FALSE,
  #adding "+" sign
  row.pos = c(1, 16:19),
  col.pos =  c(3),
  cell.text = c("+"),
  text.col = 'black',
  text.cex = 1
);

covariate_3

#############################################################
#making the legends for covariate on the left side

Q4.legends <- list(
  #legend for Patient ID
  legend = list(
    colours = c( 'blue','purple','green', 'orange', 'yellow', 'black', 'wheat4','green4','grey', 'red4') ,
    labels = c ('CPCG0001', 'CPCG0003', 'CPCG0006','CPCG0020','CPCG0042','CPCG0099','CPCG0102','CPCG0103','CPCG0183','CPCG0184'),
    title = expression(underline('Patient ID')),
    size = 3,
    continuous = FALSE,
    height = 1.5,
    pos.x = 0,
    pos.y = 0.3
  ),
  #legend for Cohort
  legend = list(
    colours = c('pink', 'royalblue'),
    labels = c('Sx', 'Bx'),
    title = expression(underline('Cohort')),
    continuous = FALSE,
    size = 3,
    height = 1.5,
    pos.x = 0,
    pos.y = 0.3
  ),
  #legend for Gleason Score
  legend = list(
    colours = c('yellow', 'orange', 'red'),
    labels = c('3+4','4+3', '4+4'),
    title = expression(underline('Gleason Score')),
    size = 3,
    continuous = FALSE,
    height = 1.5,
    pos.x = 0,
    pos.y = 0.3
  ),
  #legend for Tissue Type
  legend = list(
    colours = c(colours()[532], colours()[557]),
    labels = c('FFPE', 'Frozen'),
    title = expression(underline('Tissue Type')),
    continuous = FALSE,
    size = 3,
    height = 1.5,
    pos.x = 0,
    pos.y = 0.3
  ),
  #legend for Publication
  legend = list(
    colours = c('red','deepskyblue3', 'mediumseagreen'),
    labels = c('Baca', 'Weischenfeldt', 'Berger'),
    title = expression(underline('Publication')),
    size = 3,
    continuous = FALSE,
    height = 1.5,
    pos.x = 0,
    pos.y = 0.3
  )
);

covariate.legend.grob <- legend.grob(
  legends = Q4.legends,
  title.just = 'left',
  label.cex = 0.85,
  title.cex = 1
);

#############################################################
#making the multiplot
Q4.multiplot <-BoutrosLab.plotting.general::create.multiplot(
  plot.objects =  list(covariate_2, mainQ4_heatmap, covariate_3, covariate_1, Q4_bar_plot),
  panel.heights = c(0.3, 0.15, 1, 0.05),
  panel.widths = c(1, 0.075),
  plot.layout = c(2,4),
  layout.skip = c(FALSE, TRUE,FALSE, FALSE, FALSE, TRUE, FALSE, TRUE),
  yaxis.lab = list(NULL, c(substring(sample.order, 9, 10)), NULL, NULL, c(0, 0.5, 1)),
  yat = list(NULL, seq(1,28, 1), NULL,NULL,c(0, 0.5, 1)),
  y.spacing = c(0.2,-0.5,-0.5),
  x.spacing = c(0.2,0.1),
  yaxis.cex = 0.8,
  xaxis.cex = 0.8,
  xlab.to.xaxis.padding = 0,
  legend = list(left = list(fun = covariate.legend.grob)),
  print.new.legend = TRUE,
  yaxis.tck =  1,
  xaxis.tck = 0,
  merge.legends = TRUE,
  y.relation = 'free',
  x.relation = 'sliced',
  ylab.label = c( 'Fraction','', '', '', ''),
  ylab.cex = 1,
  resolution = 50)

Q4.multiplot
