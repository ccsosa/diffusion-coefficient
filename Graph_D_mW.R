require(ggplot2);require(reshape2);require(tidyr);require(ggpubr)

out_dir <- "E:/JAVERIANA/ASIGNATURAS/2019_2/ROTACION_DROCHSS"
  
pacman::p_load(googlesheets,dplyr)#Chemmine

#Getting Google sheets list in a personal Gmail account
(my_sheets <- gs_ls())
my_sheets %>% glimpse()

##Calling up the google sheet table desired

Sheet <- gs_title("Diffussion_coefficient_2019_09_02")
gs_ws_ls(Sheet)
#Reading table
dat <- gs_read(ss = Sheet, ws = "data")


# dat <- dat[which(!is.na(dat$`D(Cm2/s) (10-5 CM2 s-1)`)),]
dat <- dat[which(!is.na(dat$`D(Cm2/s)`)),]
dat <- dat[which(!is.na(dat$`Mm/U`)),]
# plot(dat$`Mm/U`,dat$`D(Cm2/s) (10-5 CM2 s-1)`)


# dat$Family <- NA
# xg <- split(dat$Family, 1:6)
# 
# for(i in 1:length(xg)){
#   xg[[i]] <- rep(i,length(xg[[i]]))
# }
# 
# xg <- do.call(c,xg)
# dat$Family <- xg


x <- as.data.frame.table(table(dat$Solvent))


families <- unique(dat$Family)

dat_to_graph_original <- dat[,c(1,2,4,6,7)]# 87 is solvent

unique(dat_to_graph_original$Solvent)
colnames(dat_to_graph_original) <- c("Solute","Family","Mm","D","Solvent")
dat_to_graph_original <- dat_to_graph_original[which(!is.na(dat_to_graph_original$Solvent)),]

#dat_to_graph$Family <- factor(dat_to_graph$Family,levels = 1:6)
dat_to_graph_original$Family <- factor(dat_to_graph_original$Family,levels = unique(dat_to_graph_original$Family))
dat_to_graph_original <- dat_to_graph_original[which(!is.na(dat_to_graph_original$Solvent)),]
dat_to_graph_original$Solvent <- factor(dat_to_graph_original$Solvent,levels = unique(dat_to_graph_original$Solvent))

# for(i in 1:length(families)){
#dat_to_graph <- dat_to_graph_original[which(dat_to_graph_original$Solvent=="DMSO"),]
####################################
summ_fam <- data.frame(compound= names(table(dat_to_graph_original$Family)),count= table(dat_to_graph_original$Family))
summ_fam <- summ_fam[,-2]
summ_fam <- summ_fam[which(summ_fam$count.Freq >=5),]
####################################
dat_to_graph <- subset(dat_to_graph_original,Family %in% summ_fam$compound)
#unique(dat_to_graph$Family)
####################################
x_graph <- ggscatter(dat_to_graph, x = "Mm", y = "D", size = 2,fill =  "Family",
            # title = "Solvent",
          rug = F,# Add marginal rug
          facet.by="Family",
          short.panel.labs=F,  
         # ylab = expression(paste(Dx10^{-5},"/cm"^{2},"S"^{-1})),#cm^{2}s0^{-1}
         ylab = expression(paste(D,": cm"^{2},"S"^{-1})),#cm^{2}s0^{-1} 
         xlab = expression(paste("Mm/u"," ", "g mol"^{-1})),
          conf.int=T,add = "reg.line", repel = TRUE,
          color = "black", palette = "Dark2") +
  stat_cor(aes(color = "black"), method = "spearman")+
  rremove("legend")+
  grids(linetype = "solid",color = "gray94",size = 1)
x_graph
####################################
####################################
####################################
####################################
summ_fam_or <- data.frame(compound= names(table(dat_to_graph_original$Solvent)),count= table(dat_to_graph_original$Solvent))
summ_fam_or <- summ_fam_or[,-2]
summ_fam_or <- summ_fam_or[which(summ_fam_or$count.Freq >=5),]
dat_to_graph_original <- subset(dat_to_graph_original,Solvent %in% summ_fam_or$compound)

####################################
x_graph_or <- ggscatter(dat_to_graph_original, x = "Mm", y = "D", size = 2,fill =  "Solvent",
                     rug = F,# Add marginal rug
                     facet.by="Solvent",
                     short.panel.labs=F,
                     
                     # ylab = expression(paste(Dx10^{-5},"/cm"^{2},"S"^{-1})),#cm^{2}s0^{-1}
                     ylab = expression(paste(D,": cm"^{2},"S"^{-1})),#cm^{2}s0^{-1} 
                     xlab = expression(paste("Mm/u"," ", "g mol"^{-1})),
                     conf.int=T,add = "reg.line", repel = TRUE,
                     color = "black", palette = "Dark2") +
  stat_cor(aes(color = "black"), method = "spearman")+
  rremove("legend")+
  grids(linetype = "solid",color = "gray94",size = 1)
x_graph_or

###GRAPH PER SOLVENT AND FAMILIES
#BOXPLOT

x_graph_box <- ggboxplot(dat_to_graph_original, x = "Solvent", y = "D", size = 2,fill =  "Solvent",
                        rug = F,# Add marginal rug
                        #facet.by="Solvent",
                        short.panel.labs=F,
                        
                        # ylab = expression(paste(Dx10^{-5},"/cm"^{2},"S"^{-1})),#cm^{2}s0^{-1}
                        ylab = expression(paste(D,": cm"^{2},"S"^{-1})),#cm^{2}s0^{-1} 
                        xlab = "",
                        conf.int=T,add = "reg.line", repel = TRUE,
                        color = "black", palette = "Dark2") +
  stat_cor(aes(color = "black"), method = "spearman")+
  rremove("legend")+
  grids(linetype = "solid",color = "gray94",size = 1)+
  theme(axis.text.x = element_text(angle = 90) )
x_graph_box



x <- as.data.frame.table(table(dat$Solvent))


x2 <- subset(dat,Solvent %in% summ_fam_or$compound)


solv <- unique(x2$Solvent)
xgraph_solv <- list()
for(i in 1:length(solv)){
  cat(i,"\n")
  x2i <- x2
x2a <- subset(x2i,Solvent  %in% solv[[i]])
x2a <- x2a[,c(1,2,4,6,7)]; x2a <- as.data.frame(x2a)
colnames(x2a) <- c("Solute","Family","M","D","Solvent")

x3 <- ggscatter(x2a, x = "M", y = "D", size = 2,fill =  "Family",
                        rug = F,# Add marginal rug
                        facet.by="Family",
                        short.panel.labs=F,
                        title=solv[[i]],
                        # ylab = expression(paste(Dx10^{-5},"/cm"^{2},"S"^{-1})),#cm^{2}s0^{-1}
                        ylab = expression(paste(D,": cm"^{2},"S"^{-1})),#cm^{2}s0^{-1} 
                        xlab = expression(paste("Mm/u"," ", "g mol"^{-1})),
                        conf.int=T,add = "reg.line", repel = TRUE,
                        color = "black", palette = "Dark2") +
  stat_cor(aes(color = "black"), method = "spearman")+
  rremove("legend")+
  grids(linetype = "solid",color = "gray94",size = 1)+

ggsave(paste0(out_dir,"/","graphics","/",solv[[i]],"_",Sys.Date(),".pdf"),x3,units="in",width=20,height=10,scale=2,dpi=600)
xgraph_solv[[i]] <- x3

}
  


x_author <- as.data.frame.table(table(dat$Reference))
