#require(ggplot2);require(reshape2);require(tidyr);require(ggpubr)

out_dir <- "E:/JAVERIANA/ROTACION/Graficos"

pacman::p_load(googlesheets,dplyr,reshape2,ggpubr,ggplot2,cluster,factoextra,BBmisc,readxl,openxlsx)#Chemmine

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
dat <- dat[which(!is.na(dat$species)),]
dat <- dat[which(dat$Solvent!="Hydroquinone"),]
dat <- dat[which(!is.na(dat$Molecular_density)),]

# plot(dat$`Mm/U`,dat$`D(Cm2/s) (10-5 CM2 s-1)`)
# unique(dat$Reactive)
x <- as.data.frame.table(table(dat$Solvent))
colnames(dat)
dat_to_graph_original <- dat[,c(1,2,5,6,11,
                                12,13,9)]# 87 is solvent
colnames(dat_to_graph_original) <- c("id","Solute","Classification","Mm","D",
                                     "Solvent","Viscocity","Molecular_density")
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
scaleFUN <- function(x) sprintf("%.2f", x)
scaleFUN_Y <- function(y) sprintf("%.1f", y)

y_breaks1 <- 0.2 ;y_breaks2 <- 1.5 ;

#i <- 1
sub_filtered <- lapply(1:length(Solvent_list),function(i){
  cat(paste0("i: ",i," | plot for: ",Solvent_list[[i]]),"\n")
  x_sub <- subset(dat_to_graph_original,Solvent %in% Solvent_list[[i]])
  x_sub$D <- x_sub$D*100000
  if(Solvent_list[[i]]=="Dimethyl sulfoxide"){
    y_breaks = y_breaks2
  } else {
    y_breaks = y_breaks1
  }
  
  # x_mean <- mean(x_sub$D,na.rm=T)
  # x_sd <- sd(x_sub$D,na.rm=T)
  # x_quantile <- quantile(x_sub$D)
  # # x_max <- x_mean +(2*x_sd)
  # # x_min <- x_mean -(2*x_sd)
  # x_min = x_quantile[2]
  # 
  # x_max= x_quantile[4]
  x_sub_filtered  <- x_sub
  # x_sub_filtered <- x_sub[which(x_sub$D <x_max & x_sub$D >x_min),]
  x_sub_filtered <- x_sub_filtered[which(x_sub$Mm <400 & x_sub$Mm >80),]
  x_sub_final <- x_sub_filtered
  
  
  x_mean <- mean(x_sub_filtered$Molecular_density,na.rm=T)
  x_sd <- sd(x_sub_filtered$Molecular_density,na.rm=T)
  x_quantile <- quantile(x_sub_filtered$Molecular_density)
   x_max <- x_mean +(2*x_sd)
   x_min <- x_mean -(2*x_sd)
  
   x_median <- x_quantile[3]
   x_median_min <- x_quantile[2]
   x_median_max <- x_quantile[4]
   if(nrow(x_sub_filtered) >5 ){ #nrow(x_sub_filtered)!=0 | nrow
    #plot(x_sub$Mm,x_sub$D,col="blue",pch=18)
    #points(x_sub_filtered$Mm,x_sub_filtered$D,col="red",pch=19)
    
    #plot(x_sub_final$Mm,x_sub_final$D,col=x_sub_final$cluster)
    #  x_graph <- ggscatter(x_sub, x = "Mm", y = "D",
    x_graph <- ggscatter(x_sub_filtered, x = "Mm", y = "Molecular_density",
                         
                         size = 20,
                         alpha = 0.5,
                         fill =  "black",
                         ggtheme = theme_light(),
                         title = Solvent_list[[i]],
                         rug = F,# Add marginal rug
                         #facet.by="Solvent",
                         short.panel.labs=F,
                         #dot.size = 800,
                         legend.title="",
                         cor.coef = F,
                         
                         label = "Solute",
                         font.label=c(48, "plain"),
                         
                         xlim=c(min(x_sub_filtered$Mm,na.rm=T),max(x_sub_filtered$Mm,na.rm=T)),
                         ylim=c(0,max(x_sub_filtered$Molecular_density,na.rm=T)),
                         # ylab = expression(paste(Dx10^{-5},"/cm"^{2},"S"^{-1})),#cm^{2}s0^{-1}
                         xlab =expression(paste("Molecular weight "," ", "(g mol"^{-1},")")),
                         ylab = expression(paste("Molecular density ","(\u212b"^{3},")")),
                          repel = TRUE,
                         color = "black", palette = "Dark2") +
      # stat_cor(
      #   aes(color="black",label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), 
      #   method = "pearson",size=20,
      #   label.x.npc=0.7,
      #   label.y.npc=0.95,#,#,
      #   # label.x = 3
      # )+# +
      #stat_regline_equation(aes(color=Classification))
      grids(axis = "xy",linetype = "dashed")
    x_graph  <- x_graph +
    # geom_line() +
    geom_hline(yintercept=x_mean, color='green', size=2) +
    geom_hline(yintercept=x_max, color='red', size=2) +
    geom_hline(yintercept=x_min, color='red', size=2)
    
    x_graph <- x_graph +
      bgcolor("#FFFFFF")+
      theme(panel.background = element_rect(fill = "white"),
            text=element_text(size=70),
            axis.text.x  =element_text(size=70,angle = 90),
            axis.text.y  = element_text(size=70,colour="black"),
            legend.title=element_text(size=70,colour="black"),
            axis.line = element_line(color = "black",size = 1,linetype = "solid"))+
      scale_x_continuous(breaks = seq(min(x_sub_filtered$Mm,na.rm=T), max(x_sub_filtered$Mm,na.rm=T), 20),labels=scaleFUN)+ 
      scale_y_continuous(breaks = seq(0, max(x_sub_filtered$Molecular_density,na.rm=T), 0.2),labels=scaleFUN_Y) #+
    #   theme_light()
    #+
    # xlim(c(min(x_sub$Mm),max(x_sub$Mm)))#+
    # ylim(c(min(x_sub$D),max(x_sub$D)))
    #stat_cor(aes(color = "black"), method = "pearson")#+
    #rremove("legend")+
    #grids(linetype = "solid",color = "gray94",size = 1)
    
    #x_graph
    ggsave(paste0(out_dir,"/","fix_20190930","/","MM_DENS","/mean/",as.character(Solvent_list[[i]]),"_MM_DENS_",
                  Sys.Date(),".png"),x_graph,units="in",width=22,height=19,scale=2,dpi=200)
    
    
    ##########MAD
    
    x_graph2 <- ggscatter(x_sub_filtered, x = "Mm", y = "Molecular_density",
                         
                         size = 20,
                         alpha = 0.5,
                         fill =  "black",
                         ggtheme = theme_light(),
                         title = Solvent_list[[i]],
                         rug = F,# Add marginal rug
                         #facet.by="Solvent",
                         short.panel.labs=F,
                         #dot.size = 800,
                         legend.title="",
                         cor.coef = F,
                         
                         label = "Solute",
                         font.label=c(48, "plain"),
                         
                         xlim=c(min(x_sub_filtered$Mm,na.rm=T),max(x_sub_filtered$Mm,na.rm=T)),
                         ylim=c(0,max(x_sub_filtered$Molecular_density,na.rm=T)),
                         # ylab = expression(paste(Dx10^{-5},"/cm"^{2},"S"^{-1})),#cm^{2}s0^{-1}
                         xlab =expression(paste("Molecular weight "," ", "(g mol"^{-1},")")),
                         ylab = expression(paste("Molecular density ","(\u212b"^{3},")")),
                         repel = TRUE,
                         color = "black", palette = "Dark2") +
      # stat_cor(
      #   aes(color="black",label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), 
      #   method = "pearson",size=20,
      #   label.x.npc=0.7,
      #   label.y.npc=0.95,#,#,
      #   # label.x = 3
      # )+# +
      #stat_regline_equation(aes(color=Classification))
      grids(axis = "xy",linetype = "dashed")
    x_graph2  <- x_graph2 +
      # geom_line() +
      geom_hline(yintercept=x_median, color='green', size=2) +
      geom_hline(yintercept=x_median_max, color='red', size=2) +
      geom_hline(yintercept=x_median_min, color='red', size=2)
    
    x_graph2 <- x_graph2 +
      bgcolor("#FFFFFF")+
      theme(panel.background = element_rect(fill = "white"),
            text=element_text(size=70),
            axis.text.x  =element_text(size=50,angle = 90),
            axis.text.y  = element_text(size=70,colour="black"),
            legend.title=element_text(size=70,colour="black"),
            axis.line = element_line(color = "black",size = 1,linetype = "solid"))+
      scale_x_continuous(breaks = seq(min(x_sub_filtered$Mm,na.rm=T), max(x_sub_filtered$Mm,na.rm=T), 20),labels=scaleFUN)+ 
      scale_y_continuous(breaks = seq(0, max(x_sub_filtered$Molecular_density,na.rm=T), 0.2),labels=scaleFUN_Y) #+
    #   theme_light()
    #+
    # xlim(c(min(x_sub$Mm),max(x_sub$Mm)))#+
    # ylim(c(min(x_sub$D),max(x_sub$D)))
    #stat_cor(aes(color = "black"), method = "pearson")#+
    #rremove("legend")+
    #grids(linetype = "solid",color = "gray94",size = 1)
    
    #x_graph
    ggsave(paste0(out_dir,"/","fix_20190930","/","MM_DENS","/median/",as.character(Solvent_list[[i]]),"_MM_DENS_MEDIAN_",
                  Sys.Date(),".png"),x_graph,units="in",width=22,height=19,scale=2,dpi=200)
    
    
     } else {
    cat(paste0("NO RECORDS FOR: ",Solvent_list[[i]]),"\n")
    x_sub_final <- NULL
    
  }
  return()
})
#https://www.datanovia.com/en/lessons/k-means-clustering-in-r-algorith-and-practical-examples/#targetText=K%2Dmeans%20algorithm%20requires%20users,of%20clusters%20to%20be%20produced.
