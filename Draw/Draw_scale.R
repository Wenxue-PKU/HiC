# color scaleを描画する

#==========================================
# mat lab color
#==========================================
c <- colorRampPalette(c("#00007F", "blue", "#007FFF", 
                        "cyan","#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))( 100 )

#==========================================
# gentle color
#==========================================
suppressWarnings(library("RColorBrewer"))
c <- rev(colorRampPalette(brewer.pal(10, "RdBu"))(100))  # gentle color


#==========================================
# 真ん中が白でない青から赤
#==========================================
c <- colorRampPalette(c(rgb(85,72,193, max=255), rgb(88,76,196, max=255), rgb(90,79,199, max=255), rgb(92,83,202, max=255), rgb(94,87,205, max=255), rgb(96,90,208, max=255), rgb(98,94,211, max=255), rgb(100,97,214, max=255), rgb(103,101,216, max=255), rgb(105,104,219, max=255), rgb(107,108,221, max=255), rgb(109,111,224, max=255), rgb(111,114,226, max=255), rgb(114,118,229, max=255), rgb(116,121,231, max=255), rgb(118,124,233, max=255), rgb(120,128,235, max=255), rgb(122,131,237, max=255), rgb(125,134,239, max=255), rgb(127,137,240, max=255), rgb(129,140,242, max=255), rgb(131,144,244, max=255), rgb(134,147,245, max=255), rgb(136,150,246, max=255), rgb(138,153,248, max=255), rgb(140,156,249, max=255), rgb(143,158,250, max=255), rgb(145,161,251, max=255), rgb(147,164,252, max=255), rgb(149,167,253, max=255), rgb(152,169,253, max=255), rgb(154,172,254, max=255), rgb(156,175,254, max=255), rgb(159,177,255, max=255), rgb(161,180,255, max=255), rgb(163,182,255, max=255), rgb(165,184,255, max=255), rgb(168,187,255, max=255), rgb(170,189,255, max=255), rgb(172,191,255, max=255), rgb(174,193,255, max=255), rgb(176,195,254, max=255), rgb(179,197,254, max=255), rgb(181,199,253, max=255), rgb(183,201,253, max=255), rgb(185,203,252, max=255), rgb(187,204,251, max=255), rgb(189,206,250, max=255), rgb(192,207,249, max=255), rgb(194,209,248, max=255), rgb(196,210,246, max=255), rgb(198,211,245, max=255), rgb(200,213,244, max=255), rgb(202,214,242, max=255), rgb(204,215,240, max=255), rgb(206,216,239, max=255), rgb(208,217,237, max=255), rgb(209,218,235, max=255), rgb(211,218,233, max=255), rgb(213,219,231, max=255), rgb(215,219,229, max=255), rgb(217,220,227, max=255), rgb(218,220,224, max=255), rgb(220,221,222, max=255), rgb(222,220,219, max=255), rgb(224,219,216, max=255), rgb(225,218,214, max=255), rgb(227,217,211, max=255), rgb(229,216,208, max=255), rgb(230,215,205, max=255), rgb(232,213,202, max=255), rgb(233,212,199, max=255), rgb(234,210,196, max=255), rgb(236,209,193, max=255), rgb(237,207,190, max=255), rgb(238,205,187, max=255), rgb(239,203,183, max=255), rgb(239,201,180, max=255), rgb(240,199,177, max=255), rgb(241,197,174, max=255), rgb(242,195,171, max=255), rgb(242,193,168, max=255), rgb(242,191,165, max=255), rgb(243,188,161, max=255), rgb(243,186,158, max=255), rgb(243,183,155, max=255), rgb(243,181,152, max=255), rgb(243,178,149, max=255), rgb(243,175,146, max=255), rgb(243,173,142, max=255), rgb(243,170,139, max=255), rgb(242,167,136, max=255), rgb(242,164,133, max=255), rgb(241,161,130, max=255), rgb(241,158,127, max=255), rgb(240,155,124, max=255), rgb(239,152,121, max=255), rgb(239,148,118, max=255), rgb(238,145,115, max=255), rgb(237,142,111, max=255), rgb(235,138,109, max=255), rgb(234,135,106, max=255), rgb(233,131,103, max=255), rgb(232,128,100, max=255), rgb(230,124,97, max=255), rgb(229,120,94, max=255), rgb(227,117,91, max=255), rgb(225,113,88, max=255), rgb(224,109,85, max=255), rgb(222,105,83, max=255), rgb(220,101,80, max=255), rgb(218,97,77, max=255), rgb(216,93,75, max=255), rgb(214,89,72, max=255), rgb(212,84,69, max=255), rgb(209,80,67, max=255), rgb(207,76,64, max=255), rgb(205,71,62, max=255), rgb(202,66,59, max=255), rgb(200,61,57, max=255), rgb(197,56,55, max=255), rgb(194,51,52, max=255), rgb(192,45,50, max=255), rgb(189,39,48, max=255), rgb(186,33,46, max=255), rgb(183,25,44, max=255),
                              rgb(180,15,41, max=255), rgb(177,1,39, max=255)), space="Lab")( 100 )


#==========================================
# red color
#==========================================
suppressWarnings(library("RColorBrewer"))
red <- brewer.pal(9, "Reds")[-1]
c <- colorRampPalette(red)(100) 



#==========================================
# 出力
#==========================================
data <- matrix(1:100, ncol=1)

FILE_out <- file.choose()


# 軸無しの場合
png(filename=FILE_out, width=200, height=45, units="px", bg="transparent")
par(oma=c(0,0,0,0), mar=c(0,0,0,0))
image(data, col=c, axes=F)
dummy <- dev.off()

# 軸をつける場合
png(filename=FILE_out, width=200, height=60, units="px", bg="transparent")
par(oma=c(0,0,0,0), mar=c(2.5,0.5,0,1))
image(data, col=c, axes=F)
min <- 0
max <- 100
haba <- (max - min)/4
label.name <- c(min, min + haba, min + 2*haba, min + 3*haba, max)
label.loc <- (label.name - min)/(max-min)
axis(1, line = 0.2, at=label.loc, labels=label.name, cex.axis=1.4)
dummy <- dev.off()



