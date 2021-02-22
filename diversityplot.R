library(ggplot2)
rm(list = ls())

data <- read.csv("~/Desktop/DataForDiversityFigure.csv", header = T)

data$Virus <- as.factor(data$Virus)
data$Honeybee <- as.factor(data$Honeybee)
data$Virus <- factor(data$Virus, levels = levels(data$Virus)[c(4,5,2,3,6,1)])

diversityplot <- ggplot(data)
diversityplot <- diversityplot + geom_pointrange(aes(x = Virus, y = data[,3], ymin = Lower, ymax = Upper, shape = Pool,color = Honeybee)) + theme_bw()
diversityplot <- diversityplot + ylab("Approximation to Watterson's estimator")
diversityplot <- diversityplot + scale_shape_manual(values=c(15:18))
diversityplot

tiff("~/Desktop/HostByLocNoLegend.tif", width = 297, height = 210, res = 300, units = "mm")
coeffplot
dev.off()
