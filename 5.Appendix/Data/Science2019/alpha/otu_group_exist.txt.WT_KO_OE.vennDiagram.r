
if (FALSE){
    install.packages("VennDiagram", repo="http://cran.us.r-project.org")
}
library(VennDiagram)


if ("pdf" == "png") {
    png(filename="result/alpha/otu_group_exist.txt.WT_KO_OE.vennDiagram.png", width=3, height=3,
        res=300, units="cm")
} else if ("pdf" == "eps") {
    postscript(file="result/alpha/otu_group_exist.txt.WT_KO_OE.vennDiagram.eps", onefile=FALSE, horizontal=FALSE,
               paper="special", width=10, height = 12, pointsize=10)
} else if ("pdf" == "pdf") {
    pdf(file="result/alpha/otu_group_exist.txt.WT_KO_OE.vennDiagram.pdf", width=3, height=3)
}else if ("pdf" == "svg") {
    svg(filename="result/alpha/otu_group_exist.txt.WT_KO_OE.vennDiagram.svg", width=3, height=3,
        pointsize=10)
} else {
    print("This format is currently unsupported. Please check the file <Rplots.pdf> in current directory.")
}


if (! FALSE) {

    data <- read.table(file="result/alpha/otu_group_exist.txt", sep="\t", quote="")

    num <- 0

    if("WT" != "CHENTONG"){
        WT <- data[grepl("\\<WT\\>",data[,2]),1]
        num <- num + 1
    }

    if("KO" != "CHENTONG"){
        KO <- data[grepl("\\<KO\\>",data[,2]),1]
        num <- num + 1
    }

    if("OE" != "CHENTONG"){
        OE <- data[grepl("\\<OE\\>",data[,2]),1]
        num <- num + 1
    }

    if("CHENTONG" != "CHENTONG"){
        CHENTONG <- data[grepl("\\<CHENTONG\\>",data[,2]),1]
        num <- num + 1
    }

    if("CHENTONG" != "CHENTONG"){
        CHENTONG <- data[grepl("\\<CHENTONG\\>",data[,2]),1]
        num <- num + 1
    }

    color_v <- c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3")[1:num]
    #label.col = c("orange", "white", "darkorchid4", "white", "white", "white", "white", "white", "darkblue", "white", "white", "white", "white", "darkgreen", "white"),

    if(num == 5){
        p <- venn.diagram(
            x = list(WT=WT, CHENTONG=CHENTONG,
                     CHENTONG=CHENTONG, KO=KO,
                     OE=OE),
            filename = NULL, col = "black", lwd = 1,
            fill = color_v,
            alpha = 0.50, main="",
            label.col = c("black"),
            cex = 1, fontfamily = "Helvetica",
            cat.col = color_v,cat.cex = 0.6, margin=0.3,
            cat.fontfamily = "Helvetica"
        )
    }else if(num == 4){
        p <- venn.diagram(
            x = list(WT=WT, CHENTONG=CHENTONG, KO=KO,
                     OE=OE),
            filename = NULL, col = "black", lwd = 1,
            fill = color_v,
            alpha = 0.50, main="",
            label.col = c("black"),
            cex = 1, fontfamily = "Helvetica",
            cat.col = color_v, cat.cex = 0.8, margin=0.2,
            cat.fontfamily = "Helvetica",
        )
    } else if (num==3) {
        p <- venn.diagram(
            x = list(WT=WT, KO=KO, OE=OE),
            filename = NULL, col = "transparent",
            fill = color_v,
            alpha = 0.50, main="",
            label.col = c("black", "black", "black", "black", "black", "black", "black"),
            cex = 1, fontfamily = "Helvetica", cat.default.pos="text",
            cat.pos=0,  magrin=0.1,
            cat.col = color_v,cat.cex = 1,cat.fontfamily = "Helvetica"
        )
    } else if (num==2) {
        p <- venn.diagram(
            x = list(WT=WT, KO=KO),
            filename = NULL, col = "transparent",
            fill = color_v,
            alpha = 0.50, main="",
            label.col = c("black"),
            cex = 1, fontfamily = "Helvetica",
            cat.default.pos="outer",
            cat.pos=0, margin=0.1,
            cat.col = color_v,cat.cex = 1,cat.fontfamily = "Helvetica"
        )
    }
    grid.draw(p)
} else {
    #---venn plot for given numbers---------
    numList <- c()
    labelList <- c()

    num <- length(labelList)
    color_v <- c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3")[1:num]

    if (num==2) {
        draw.pairwise.venn(area1=numList[1], area2=numList[2],
                           cross.area=numList[3], category=labelList, lwd=rep(1,1),
                           lty=rep(2,2), col="transparent", fill=color_v,
                           cat.col=color_v)
    } else if (num==3) {
        draw.triple.venn(area1=numList[1], area2=numList[2],
                         area3=numList[3], n12=numList[4], n23=numList[5],
                         n13=numList[6], n123=numList[7],
                         category=labelList, col="transparent", fill=color_v,
                         cat.col=color_v, reverse=FALSE)
    }else if (num==4) {
        draw.quad.venn(area1=numList[1], area2=numList[2],
                       area3=numList[3], area4=numList[4], n12=numList[5],
                       n13=numList[6], n14=numList[7], n23=numList[8],
                       n24=numList[9], n34=numList[10], n123=numList[11],
                       n124=numList[12], n134=numList[13], n234=numList[14],
                       n1234=numList[15],
                       category=labelList, col="transparent", fill=color_v,
                       cat.col=color_v, reverse=FALSE)
    }else if (num==5){
        #draw.quintuple.venn()
    }
}

dev.off()

