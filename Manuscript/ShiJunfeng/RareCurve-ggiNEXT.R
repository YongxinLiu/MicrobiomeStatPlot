if (!requireNamespace("iNEXT", quietly = TRUE))
  install.packages("iNEXT")
if (!requireNamespace("ggplot2", quietly = TRUE))
  install.packages("ggplot2")
if (!requireNamespace("tidyverse", quietly = TRUE))
  install.packages("tidyverse")
library(iNEXT)
library(ggplot2)
suppressWarnings(suppressMessages(library(tidyverse)))
d <- read.table("otu_table.txt", sep="\t", dec=".", header=T)
d<-subset(d, select = -X )
out <- iNEXT(sapply(d, as.numeric), q=0, datatype="abundance", knots=100)
ggiNEXT(out, type=1, se = TRUE, grey =FALSE)+
  ylim(c(0,750)) +
  xlim(c(0,90000)) + scale_color_manual(breaks = c("a1", "a2", "a3","b1","b2","b3"),
                                        values=c("#C77CFF","#C77CFF","#C77CFF","#00BA38","#00BA38","#00BA38"))