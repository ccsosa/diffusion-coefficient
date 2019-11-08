#require(ggplot2);require(reshape2);require(tidyr);require(ggpubr)

out_dir <- "E:/JAVERIANA/ROTACION/Graficos"

pacman::p_load(googlesheets,dplyr,reshape2,ggpubr,ggplot2,cluster,factoextra,BBmisc,readxl,openxlsx)#Chemmine


inDir <- "E:/JAVERIANA/ROTACION"

data <- read.xlsx(paste0(inDir,"/","sub_filtered_D_CALC.xlsx"),sheet = 1)

data$D <- data$D*100000
data$D_CALC <- data$D_CALC*100000
plot(data$D,data$D_CALC)

data2 <- data[which(data$Solvent!="Dimethyl sulfoxide"),]
x_graph_RMSE <- ggscatter(data2, x = "D", y = "D_CALC",
                     #colour="color_1",
                     size = 5,
                     alpha = 0.5,
                     fill =  "Solvent",
                     ggtheme = theme_light(),
                     title = "",
                     rug = F,# Add marginal rug
                     #facet.by="Solvent",
                     short.panel.labs=F,
                     #dot.size = 800,
                     legend.title="",
                     cor.coef = F,
                     
                     #label = "Solute",
                     #font.label=c(48, "plain"),
                     
                     # xlim=c(min(x_sub_filtered$Mm,na.rm=T)-10,max(x_sub_filtered$Mm,na.rm=T)+10),
                     # ylim=c(min(x_sub_filtered$D,na.rm=T),max(x_sub_filtered$D,na.rm=T)),
                     # # ylab = expression(paste(Dx10^{-5},"/cm"^{2},"S"^{-1})),#cm^{2}s0^{-1}
                     ylab = expression(paste("observed D (cm"^{2},"S"^{-1},") x 10"^{-5})),#cm^{2}s0^{-1} 
                     xlab = expression(paste("predicted D (cm"^{2},"S"^{-1},") x 10"^{-5})),#cm^{2}s0^{-1} 
                     #xlab = expression(paste("Molecular weight "," ", "(g mol"^{-1},")")),
                     conf.int=T,add = "reg.line", repel = TRUE,
                     color = "Solvent",
                     palette = "Paired")+
  theme(panel.background = element_rect(fill = "white"),
        text=element_text(size=40),
        axis.text.x  =element_text(size=50,angle = 90),
        axis.text.y  = element_text(size=70,colour="black")
        #legend.title=element_text(size=70,colour="black")
        )
  
  
  # + #"Dark2"
# 
#   stat_cor(
#     aes(color=Solvent,label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
#     method = "pearson",size=10,
#     label.x.npc=0.,
#     label.y.npc=0.95,#,#,
#     # label.x = 3
#   )+# +
#   #stat_regline_equation(aes(color=Classification))
#   grids(axis = "xy",linetype = "dashed")

x_graph_RMSE

x_graph_RMSE2 <- ggscatter(data, x = "D", y = "D_CALC",
                          #colour="color_1",
                          size = 5,
                          alpha = 0.5,
                          fill =  "Solvent",
                          ggtheme = theme_light(),
                          title = "",
                          rug = F,# Add marginal rug
                          #facet.by="Solvent",
                          short.panel.labs=F,
                          #dot.size = 800,
                          legend.title="",
                          cor.coef = F,
                          
                          #label = "Solute",
                          #font.label=c(48, "plain"),
                          
                          # xlim=c(min(x_sub_filtered$Mm,na.rm=T)-10,max(x_sub_filtered$Mm,na.rm=T)+10),
                          # ylim=c(min(x_sub_filtered$D,na.rm=T),max(x_sub_filtered$D,na.rm=T)),
                          # # ylab = expression(paste(Dx10^{-5},"/cm"^{2},"S"^{-1})),#cm^{2}s0^{-1}
                          ylab = expression(paste("observed D (cm"^{2},"S"^{-1},") x 10"^{-5})),#cm^{2}s0^{-1} 
                          xlab = expression(paste("predicted D (cm"^{2},"S"^{-1},") x 10"^{-5})),#cm^{2}s0^{-1} 
                          #xlab = expression(paste("Molecular weight "," ", "(g mol"^{-1},")")),
                          conf.int=T,add = "reg.line", repel = TRUE,
                          color = "Solvent",
                          palette = "Paired")+
  theme(panel.background = element_rect(fill = "white"),
        text=element_text(size=40),
        axis.text.x  =element_text(size=50,angle = 90),
        axis.text.y  = element_text(size=70,colour="black")
        #legend.title=element_text(size=70,colour="black")
  )
#   stat_cor(
#     aes(color=Solvent,label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
#     method = "pearson",size=10,
#     label.x.npc=0.,
#     label.y.npc=0.95,#,#,
#     # label.x = 3
#   )+# +
#   #stat_regline_equation(aes(color=Classification))
#   grids(axis = "xy",linetype = "dashed")

x_graph_RMSE2

