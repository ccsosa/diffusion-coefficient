#require(ggplot2);require(reshape2);require(tidyr);require(ggpubr)

out_dir <- "E:/JAVERIANA/ROTACION/Graficos"

pacman::p_load(googlesheets,dplyr,reshape2,ggpubr,ggplot2)#Chemmine

#Getting Google sheets list in a personal Gmail account
(my_sheets <- gs_ls())
my_sheets %>% glimpse()
Sheet <- gs_title("Diffussion_coefficient_2019_09_02")
gs_ws_ls(Sheet)
#Reading table
dat <- gs_read(ss = Sheet, ws = "data")


# dat <- dat[which(!is.na(dat$`D(Cm2/s) (10-5 CM2 s-1)`)),]
dat <- dat[which(!is.na(dat$`D(Cm2/s)`)),]
dat <- dat[which(!is.na(dat$`Mm/U`)),]
dat <- dat[which(!is.na(dat$Reactive)),]
# plot(dat$`Mm/U`,dat$`D(Cm2/s) (10-5 CM2 s-1)`)
unique(dat$Reactive)
x <- as.data.frame.table(table(dat$Solvent))

dat_to_graph_original <- dat[,c(2,5,6,8,9)]# 87 is solvent
colnames(dat_to_graph_original) <- c("Solute","Classification","Mm","D","Solvent")
dat_to_graph_original <- dat_to_graph_original[which(!is.na(dat_to_graph_original$Solvent)),]

#dat_to_graph$Family <- factor(dat_to_graph$Family,levels = 1:6)
dat_to_graph_original$Classification <- factor(dat_to_graph_original$Classification,levels = unique(dat_to_graph_original$Classification))
dat_to_graph_original <- dat_to_graph_original[which(!is.na(dat_to_graph_original$Solvent)),]
dat_to_graph_original$Solvent <- factor(dat_to_graph_original$Solvent,levels = unique(dat_to_graph_original$Solvent))



####################################
summ_fam <- data.frame(compound= names(table(dat_to_graph_original$Classification)),count= table(dat_to_graph_original$Classification))
summ_fam <- summ_fam[,-2]
summ_fam <- summ_fam[which(summ_fam$count.Freq >=5),]
####################################
dat_to_graph <- subset(dat_to_graph_original,Classification %in% summ_fam$compound )
Solvent_list <- unique(dat_to_graph$Solvent)

#Our transformation function
scaleFUN <- function(x) sprintf("%.0f", x)
scaleFUN_Y <- function(y) sprintf("%.1f", y)

y_breaks1 <- 0.3 ;y_breaks2 <- 1.5 ;
lapply(1:length(Solvent_list),function(i){
cat(paste0("plot for:",Solvent_list[[i]]),"\n")
x_sub <- subset(dat_to_graph_original,Solvent %in% Solvent_list[[i]])
x_sub$D <- x_sub$D*100000
if(Solvent_list[[i]]=="Dimethyl sulfoxide"){
  y_breaks = y_breaks2 
} else {
  y_breaks = y_breaks1 
}
x_graph <- ggscatter(x_sub, x = "Mm", y = "D",
                     size = 20,
                     alpha = 0.5,
                     fill =  "Classification",
                      title = Solvent_list[[i]],
                     rug = F,# Add marginal rug
                     #facet.by="Solvent",
                     short.panel.labs=F,
                     #dot.size = 800,
                     legend.title="",
                     cor.coef = F,
                     
                     label = "Solute",
                     font.label=c(48, "plain"),
                     
                     xlim=c(min(x_sub$Mm)-10,max(x_sub$Mm)+10),
                     ylim=c(min(x_sub$D),max(x_sub$D)),
                     # ylab = expression(paste(Dx10^{-5},"/cm"^{2},"S"^{-1})),#cm^{2}s0^{-1}
                     ylab = expression(paste("Diffusion coefficient (cm"^{2},"S"^{-1},") x 10"^{-5})),#cm^{2}s0^{-1} 
                     xlab = expression(paste("Molecular weight "," ", "(g mol"^{-1},")")),
                     conf.int=T,add = "reg.line", repel = TRUE,
                     color = "Classification", palette = "Dark2") +
  stat_cor(
    aes(color=Classification,label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), 
        method = "pearson",size=20,
    label.x.npc=0.7,
    label.y.npc=0.95,#,#,
   # label.x = 3
  )+# +
  #stat_regline_equation(aes(color=Classification))
  grids(axis = "xy",linetype = "dashed")
  
  x_graph <- x_graph +
  bgcolor("#FFFFFF")+
  theme(panel.background = element_rect(fill = "gray90"),
        text=element_text(size=70),
        axis.text.x  =element_text(size=50,angle = 90),
        axis.text.y  = element_text(size=70,colour="black"),
        legend.title=element_text(size=70,colour="black"))+
  scale_x_continuous(breaks = seq(min(x_sub$Mm), max(x_sub$Mm)+40, 20),labels=scaleFUN)+ 
  scale_y_continuous(breaks = seq(min(x_sub$D), max(x_sub$D), y_breaks),labels=scaleFUN_Y)

 #+
# xlim(c(min(x_sub$Mm),max(x_sub$Mm)))#+
  # ylim(c(min(x_sub$D),max(x_sub$D)))
  #stat_cor(aes(color = "black"), method = "pearson")#+
  #rremove("legend")+
  #grids(linetype = "solid",color = "gray94",size = 1)

#x_graph
ggsave(paste0(out_dir,"/",as.character(Solvent_list[[i]]),"_",Sys.Date(),".png"),x_graph,units="in",width=22,height=19,scale=2,dpi=300)

})
####################################
dat_to_graph_original2 <- dat_to_graph_original
dat_to_graph_original2$Mm <- dat_to_graph_original2$Mm/100000


x_graph_box <- ggboxplot(dat_to_graph_original2, x = "Solvent", y = "D", size = 2,fill =  "Classification",
                         rug = F,# Add marginal rug
                         #facet.by="Solvent",
                         short.panel.labs=F,
                         
                         # ylab = expression(paste(Dx10^{-5},"/cm"^{2},"S"^{-1})),#cm^{2}s0^{-1}
                         ylab = expression(paste("Diffusion coefficient (cm"^{2},"S"^{-1},") x 10"^{-5})),#cm^{2}s0^{-1} 
                         xlab = "",
                         conf.int=T,add = "reg.line", repel = TRUE,
                         #color = "black", 
                         palette = "Dark2") +
  stat_cor(aes(color = "black"), method = "spearman")+
  # rremove("legend")+
  grids(linetype = "solid",color = "gray94",size = 1)+
  theme(axis.text.x = element_text(angle = 90) )

x_graph_box
