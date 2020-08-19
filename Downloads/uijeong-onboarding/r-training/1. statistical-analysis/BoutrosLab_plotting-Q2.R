#SIMPLE PLOTS: BoutrosLab.plotting.general
#=======
#  2a. Create a simple scatterplot using BoutrosLab.plotting.general.
#-Use the 'cars' dataset provided the R Datasets package
#-Set a title and x- and y-axis labels that make sense
#-Ensure that the axes ranges make sense
#-Fiddle with other parameters as needed to make the figure aesthetically pleasing (e.g. font size, line width, etc)

library(BoutrosLab.plotting.general);
library(dplyr);

head(cars)
head(mtcars)
scatterplot.mpg.displacement <- BoutrosLab.plotting.general:::create.scatterplot(
  formula = disp ~ mpg, 
  data = mtcars, main =  'Motor Trends Car Data',  
  xlab.label = 'MPG', 
  ylab.label = 'Displacement', 
  main.cex = 2, 
  xlab.cex = 1.5, 
  ylab.cex = 1.5, 
  xlimits = c(8, 38), 
  ylimits = c(50, 500), 
  xgrid.at = seq(10, 35, 5),
  col = 'blue',
  xat =  seq(10, 35, 5),
  xaxis.cex = 1, 
  yaxis.cex = 1, 
  x.spacing = 3, 
  y.spacing =3,
  resolution = 50,
  type = c('p', 'g')
);

scatterplot.mpg.displacement

#2b. Try another function now. Create a heatmap displaying data found in the 'Loblolly' dataset found in the R Datasets package. 
#-Note: you may need to reorganize the data
#-Do all formatting necessary to improve presentation (e.g. colour scheme, etc)
head(Loblolly)
Loblolly.wide <- reshape(Loblolly, idvar = 'age', timevar = 'Seed', direction = 'wide');
Loblolly.wide

BoutrosLab.plotting.general:::create.heatmap(
  x = Loblolly.wide[,-1],
  clustering.method = 'none',
  main = '  Loblolly growth from independent seeds',
  main.cex = 2,
  xlab.label = 'Time(years)',
  xlab.cex = 1.5,
  xaxis.cex = 0.75,
  ylab.label = 'Seed',
  ylab.cex = 1.5,
  yaxis.cex = 0.75,
  yaxis.lab = (c(substring(colnames(Loblolly.wide[,-1]), 8, 10))),
  xaxis.lab = c(row.names(Loblolly.wide)),
  xaxis.rot = 0,
  colour.scheme = c('white', 'firebrick'),
  colour.alpha = 'automatic',
  print.colour.key = FALSE,
  xaxis.tck = 0,
  top.padding = 2,
  bottom.padding = 2,
  left.padding = 2,
  right.padding = 2,
  resolution = 50)

#2c. Take a look at the 'ChickWeight' dataset in the R Datasets package
#-Think about which plot-type should be used to effectively highlight the relationship between the growth of chicks and their diet.
#-Create the plot using BoutrosLab.plotting.general functions
#-Note that you have many options here, including plotting all the data points, aggregating the data, controlling for certain factors, etc.

head(ChickWeight)
tail(ChickWeight)
str(ChickWeight)
names(ChickWeight)

ChickWeight
preDiet.mean <-ChickWeight %>% filter(Time == 0) %>% group_by(Diet) %>% summarize(mean=as.numeric(mean(weight, na.rm=TRUE)));
postDiet.mean <- ChickWeight %>% filter(Time == 20) %>% group_by(Diet) %>% summarize(mean=mean(weight, na.rm=TRUE));
preDiet.median <- ChickWeight %>% filter(Time == 0) %>% group_by(Diet) %>% summarize(median=median(weight, na.rm=TRUE));
postDiet.median <- ChickWeight %>% filter(Time == 20) %>% group_by(Diet) %>% summarize(median=median(weight, na.rm=TRUE));


data.merge <- rbind(preDiet.mean, postDiet.mean);
data.merge <- cbind(data.merge,c(rep(0,4), rep(20,4)));
colnames(data.merge)<- c("Diet","MeanOfWeight","Time");
data.merge

BoutrosLab.plotting.general::create.barplot(
  formula = MeanOfWeight ~ Diet, 
  data = data.merge,
  main = 'Effect of diet on chick weight',
  ylab.label = 'The mean of weight gained',
  xlab.label = 'Diet',
  ylimits = c(0,350),
  main.cex = 1.75,
  ylab.cex = 1.25,
  xlab.cex = 1.25,
  xaxis.fontface = 1,
  yaxis.fontface = 1,
  # Setting groups
  groups = Time,
  col = default.colours(12, is.greyscale = FALSE)[11:12],
  legend = list(
    inside = list(
      fun = draw.key,
      args = list(
        key = list(
          points = list(
            col = 'black',
            pch = 22,
            cex = 1.5,
            fill = default.colours(12, is.greyscale = FALSE)[11:12]
          ),  
          text = list(
            lab = c('Pre-Diet(0 days)','Post-Diet(after 20 days')
          ),
          padding.text = 2,
          cex = 0.8
        )
      ),
      x = 0.65,
      y = 0.96
    )
  ),
  description = 'Barplot created by BoutrosLab.plotting.general',
  resolution = 50
)

