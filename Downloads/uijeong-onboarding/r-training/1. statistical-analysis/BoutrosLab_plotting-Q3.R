#COMPLEX PLOTS: BoutrosLab.plotting.general
#=======
#  3. 
#- Take a look at "Q3_SampleOutput.tiff"
#- This is the figure you will be re-creating using BoutrosLab.plotting.general functions. 
#- The data used to generate the figure is found in the "Q3_SeqControl_data.tsv" file

#======= BACKGROUND INFO =======
#  The first thing you need to do is understand the figure. This is a figure used to describe the performance of SeqControl. 

#SeqControl is a framework used to make predictions in sequencing data quality. For example, a researcher may want to conduct an experiment which requires sequencing data. The researcher may begin by generating a subset of the data, and measuring metrics of the data. These metrics can be combined into a statistical model to predict features of the experiment, in order to improve the choices made when sequencing the rest of the data required. 

#An overview figure "SeqControl_Overview.tff" can be found in the outputs folder 

#The statistical model referred to in this figure is called a random forest. Essentially, this model makes a prediction regarding the quality of sequencing data that will be generated. The fraction of 'yes' vote indicates the confidence that six lanes of sequencing will achieve 50x coverage. 

#The actual observed coverage of the data is shown by the colour of the bar, with black indicating high coverage, and grey indicating lower than 50x coverage. 

#There are also five covariate bars beneath the barplot.

#The first covariate bar indicates the tumour sample from which the sequencing data came from. The second covariate bar indicates how the tumour sample was prepared - either FFPE or frozen. The remaining three covariates indicate some metric related to the data:
#  -% Bases > 0
#-Unique start points
#-Average reads/start

#======= UNDERSTANDING THE DATA =======
#  The data for this plot has been organized into a tab-delimited file. The header of the file indicates what kind of data is stored in each column.

#======= PLOTTING =======
#Overview: The main plot in this figure is made using create.barplot. The covariates are made using create.heatmap. The legend is created using legend.grob. All of these elements are combined together using create.multiplot.
library(BoutrosLab.plotting.general);
library(dplyr);

#############################################################
#Step 1. You will need to read in the data and do some reformatting. If you look at the final plot, the "yes votes" are ordered in decreasing order. 
#This is not how they are ordered in the data file. You will need to reorder the data after you read it in (don't make changes to the original data file). 
SeqControl <- read.delim("Q3/Q3_SeqControl_data.tsv", header = TRUE, sep ="\t");
head(SeqControl)

Q3_SeqControl_data <- arrange(SeqControl, desc(yes.votes));
head(Q3_SeqControl_data)

#############################################################
#Step 2. Create a heatmap displaying the cpcgene sample names. 
#The colour scheme to use for this is given in the "Q3_PlottingInfo.R" file. 
#(Note that create.multiplot can be used to adjust plot proportions at the end).
sample.colour.scheme <- c('rosybrown1', 'rosybrown4', 'red', 'darkred', 'darkorange', 'gold', 'darkolivegreen3', 'darkgreen', 'aquamarine', 'cyan4', 'dodgerblue', 'darkblue');
sample.names <- c('CPCG0003P', 'CPCG0005P', 'CPCG0007P', 'CPCG0040P', 'CPCG0047P', 'CPCG0603P', 'CPCG0098P', 'CPCG0102P', 'CPCG0103P', 'CPCG0123P', 'CPCG0183P', 'CPCG0184P');

#re-generating the cpcgene sample names.
head(Q3_SeqControl_data)
tmp.cpcg <- matrix(unlist(strsplit(as.character(Q3_SeqControl_data$CPCG), split="CPCG")));
Q3_SeqControl_cpcg <- tmp.cpcg[!(tmp.cpcg[,1] == ""),];
Q3_SeqControl_cpcg <- as.numeric(matrix(unlist(strsplit(as.character(Q3_SeqControl_cpcg), split="P"))));
Q3_SeqControl_cpcg <- as.data.frame(Q3_SeqControl_cpcg);
Q3_SeqControl_cpcg

cpcgene_neames <- create.heatmap(
  x = t(Q3_SeqControl_cpcg),
  clustering.method = 'none',
  scale.data = FALSE,
  force.grid.col = TRUE,
  colour.scheme = sample.colour.scheme,
  grid.col = TRUE,
  print.colour.key = FALSE,
  yaxis.tck = 0,
)

cpcgene_neames

#############################################################
#Step 3. Create heatmaps for each of the three metrics. 
#If you read the header names, "Average.reads.start", "Unique.start.points", and "X..Bases...0.quality" are the columns you are looking for. 
#The colour schemes for these range from white to either "deeppink", "darkblue", or "darkorange" (these are R colour names). 
#Note that you will have to show the metric values on a continuous scale. 
head(Q3_SeqControl_data)

#creating the "Average.reads.start" covariate bar
avg_reads_starts.mtx <- as.matrix(as.numeric(Q3_SeqControl_data$Average.reads.start));

Avg_reads_starts_covariate <- create.heatmap(
  x = t(avg_reads_starts.mtx),
  clustering.method = 'none',
  scale.data = FALSE,
  force.grid.col = TRUE,
  colour.scheme = c('white','deeppink'),
  grid.col = TRUE,
  print.colour.key = FALSE,
  yaxis.tck = 0,  #remove y axis ticks
)

Avg_reads_starts_covariate

#############################################################
#creating the "Unique.start.points" covariate bar
Unique.start.points.mtx <- as.matrix(as.numeric(Q3_SeqControl_data$Unique.start.points));

Unique.start.points_covariate <- create.heatmap(
  x = t(Unique.start.points.mtx),
  clustering.method = 'none',
  scale.data = FALSE,
  force.grid.col = TRUE,
  colour.scheme = c('white','darkblue'),
  grid.col = TRUE,
  print.colour.key = FALSE,
  yaxis.tck = 0, 
)

Unique.start.points_covariate

#############################################################
#creating the "X..Bases...0.quality" covariate bar
Xbases0.matrix <- as.matrix(as.numeric(Q3_SeqControl_data$X..Bases...0x));

X.Bases.0_covariate <- create.heatmap(
  x = t(Xbases0.matrix),
  clustering.method = 'none',
  scale.data = FALSE,
  force.grid.col = TRUE,
  colour.scheme = c('white','darkorange'),
  grid.col = TRUE,
  print.colour.key = FALSE,
  yaxis.tck = 0,
)

X.Bases.0_covariate

#############################################################
#Step 4. For the "FFPE" covariate bar, you'll need to know that the samples which are FFPE are samples "CPCG0102P" and "CPCG0103P". Also, the colour scheme consists of "white" and "darkslategrey".
tail(Q3_SeqControl_data)
FFPE.matrix <- as.numeric(Q3_SeqControl_data$CPCG %in% c('CPCG0102P', 'CPCG0103P'));

#creating the heatmap
FFPE_covariate <- create.heatmap(
  x = t(FFPE.matrix),
  clustering.method = 'none',
  scale.data = FALSE,
  force.grid.col = TRUE,
  colour.scheme = c('white' , 'darkslategrey'),
  total.col = 12,
  grid.col = TRUE,
  print.colour.key = FALSE,
  yaxis.tck = 0,
)

FFPE_covariate

#############################################################
#Step 5. Plot the barplot. The colour scheme consists of "black" and "grey". 
#The "outcome" column in the data file tells you what the observed coverage was for a sample. An outcome of 0 indicates coverage <50x, and 1 indicates coverage >=50x.

#splicing data frame for the barplot
Q3_bar_plot_data.df <- data.frame( 
  yes.votes = Q3_SeqControl_data$yes.votes,
  outcome = Q3_SeqControl_data$outcome,
  num = c(1:72)
);
Q3_bar_plot_data.df

#colors
Q3_bar_plot.colors <- replace(as.vector(Q3_bar_plot_data.df$outcome), which(Q3_bar_plot_data.df$outcome == 0), 'grey');
Q3_bar_plot.colors <- replace(Q3_bar_plot.colors, which(Q3_bar_plot_data.df$outcome == 1), 'black');
Q3_bar_plot.colors

# making the barplot
Q3_bar_plot <- create.barplot(
  formula = yes.votes ~ num,
  data = Q3_bar_plot_data.df,
  ylab.label = "Fraction of yes votes",
  ylimits = c(0,1.05),
  xaxis.tck = 0,
  xaxis.cex = 0,
  yaxis.cex = 1,
  xlab.cex = 0,
  ylab.cex = 1.5,
  yat = seq(0.0, 1.0, 0.25),
  sample.order = 'decreasing',
  col = Q3_bar_plot.colors, 
  abline.h = 0.50,
  abline.lwd = 1,
  abline.col = 'grey',
  abline.lty = 'dashed',
  legend = list(
    inside = list(
      fun = draw.key,
      args = list(
        key = list(
          title = expression(underline('Observed')),
          points = list(
            col = 'black',
            pch = 22,
            cex = 2,
            fill = c('grey', 'black')),
          text = list(
            lab = c('< 50X', '>= 50X')),
          cex = 1)),
      x = 0.85,  #position
      y = 0.95
      )
    ),
  resolution = 50
  );

Q3_bar_plot       

#############################################################
#Step 6. Create a legend for each of the covariates. Ensure that the legend displays either continuous or discrete schemes as appropriate. 
#The individual legends can be combined using the legend.grob function. The labels for "unique start points" can be formatted using the scientific.notation function.

#making the legends for the covariates on the right side
Q3.legends <- list(
    #legend for sample
    legend = list(
      title = expression(underline('Sample')),
      colours = sample.colour.scheme,
      labels = sample.names,
      size = 1.5,
      border = 'black'
    ),
    #legend for Sample Preparation
    legend = list(
      title = expression(underline('Sample Preparation')),
      colours = c('white','darkslategrey'),
      labels = c('Frozen', 'FFPE'),
      size = 1.5,
      border = 'black'
    ),
    #legend for % Bases > 0 quality
    legend = list(
      title = expression(underline('% Bases > 0 quality')),
      colours = c('darkorange', 'white'),
      labels = c('97.0', '83.0'),
      border = 'black',
      continuous = TRUE
    ),
    #legend for Unique start points
    legend = list(
      title = expression(underline('Unique start points')),
      colours = c('darkblue', 'white'),
      labels = scientific.notation(c(360000000, 100000000)),
      border = 'black',
      continuous = TRUE
    ),
    #legend for Average reads/start
    legend = list(
      title = expression(underline('Average reads/start')),
      colours = c('deeppink', 'white'),
      labels = c('1.190', '1.070'),
      border = 'black',
      continuous = TRUE
    )
);

covariate.legend.grob <- legend.grob(
  legends = Q3.legends,
  title.just = 'left',
  label.cex = 0.85,
  title.cex = 1
);

#############################################################
#Step 7. Create a legend for the barplot. This should be created separately for ease of placement in the final figure.
#recreating code for genenames with color codes
tmp.cpcg <- matrix(unlist(strsplit(as.character(Q3_SeqControl_data$CPCG), split="CPCG")));
Q3_SeqControl_cpcg <- tmp.cpcg[!(tmp.cpcg[,1] == ""),];
Q3_SeqControl_cpcg <- as.numeric(matrix(unlist(strsplit(as.character(Q3_SeqControl_cpcg), split="P"))));
Q3_SeqControl_cpcg <- as.data.frame(Q3_SeqControl_cpcg);
Q3_SeqControl_cpcg

#making the heatmap for the sample covariate
CPCG_covariate <- create.heatmap(
  x = t(Q3_SeqControl_cpcg),
  clustering.method = 'none',
  scale.data = FALSE,
  force.grid.col = TRUE,
  colour.scheme = sample.colour.scheme,
  grid.col = TRUE,
  print.colour.key = FALSE,
  yaxis.tck = 0,
  xaxis.tck = 0,
  legend.side = 'right'
)

CPCG_covariate

#############################################################                                                                                                                                                                                                                                                               
#Step 8. Combine all of the plots and legends together using create.multiplot.

#making the multiplot          
Q3.multiplot <- create.multiplot(
  plot.objects =  list(Avg_reads_starts_covariate,Unique.start.points_covariate,X.Bases.0_covariate, FFPE_covariate, CPCG_covariate, Q3_bar_plot),
  panel.heights = c(15, 1, 1, 1, 1, 1),
  xaxis.lab = NULL,
  yaxis.lab = list(NULL, NULL, NULL, NULL, NULL, seq(0.00, 1.00, 0.25)),
  ylab.label = c("Fraction of var votes",''),
  y.spacing = c(0.3, 0.3, 0.3, 0.3, 0.3),
  legend = list(right = list(fun = covariate.legend.grob)),
  print.new.legend = TRUE,
  xaxis.alternating = 0,
  ylab.cex = 1.5,
  xlab.to.xaxis.padding = 0,
  xaxis.tck = 0,
  yaxis.tck = c(0, 2),
  merge.legends = TRUE,
  resolution = 50
)

Q3.multiplot
