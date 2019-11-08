  #require(ggplot2);require(reshape2);require(tidyr);require(ggpubr)
  
  out_dir <- "E:/JAVERIANA/ROTACION/Graficos"
  
  pacman::p_load(googlesheets,dplyr,reshape2,ggpubr,ggplot2,plm,cluster,factoextra,BBmisc,readxl,openxlsx)#Chemmine

  
  
  
  #################################3
  colors <- data.frame( classification= c(
  "Hydroxy group",
  "Halogen",
  "Like Benzene",
  "Like Phenyl",
  "No reactive",
  "Others"
                       ), color = c( "#E51A82",
                                     "#D75700",
                                     "#8B5A00",
                                     "#7570B3",
                                     "#159B73",
                                     "#978E72"
                                     
    
              )
  )
  
  
  
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
  scaleFUN <- function(x) sprintf("%.0f", x)
  scaleFUN_Y <- function(y) sprintf("%.1f", y)
  
  y_breaks1 <- 0.2 ;y_breaks2 <- 1.5 ;
  
  #i <- 1
  sub_filtered <- lapply(1:length(Solvent_list),function(i){
    cat(paste0("i: ",i," | plot for: ",Solvent_list[[i]]),"\n")
    x_sub <- subset(dat_to_graph_original,Solvent %in% Solvent_list[[i]])
    x_sub1 <- x_sub
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
    x_sub_filtered <- as.data.frame(x_sub_filtered)
    x_sub_filtered$Classification <- factor(x_sub_filtered$Classification,levels = unique(x_sub_filtered$Classification))
    x_sub_filtered$Solvent <- factor(x_sub_filtered$Solvent,levels = unique(x_sub_filtered$Solvent))

        
    if(nrow(x_sub_filtered) >=5 ){ #nrow(x_sub_filtered)!=0 | nrow
    #plot(x_sub$Mm,x_sub$D,col="blue",pch=18)
    #points(x_sub_filtered$Mm,x_sub_filtered$D,col="red",pch=19)
    x_sub_filtered_scaled <- as.data.frame(scale(cbind(x_sub_filtered$Mm,x_sub_filtered$D)))
    x_sub_filtered_scaled$id <- x_sub_filtered$id
    x_sub_filtered_scaled <- x_sub_filtered_scaled[which(complete.cases(x_sub_filtered_scaled)),]
    #set.seed(123)
    
    n_clust <- fviz_nbclust(x_sub_filtered_scaled[,1:2],FUNcluster = kmeans,k.max = 5,nboot=100) 
    n_clust <- n_clust$data[which.max(n_clust$data$y),][1]
    
    
    km.res <- kmeans(x_sub_filtered_scaled[,1:2], as.numeric(n_clust),
                     nstart = 4,iter.max = 1000,)
    dd <- as.data.frame(cbind(x_sub_filtered_scaled, cluster = km.res$cluster))
    dd$id <- factor(as.character(dd$id))
      # plot(dd$V1,dd$V2,col=dd$cluster)
    
    x_sub_filtered$id <- factor(as.character(x_sub_filtered$id))
    x_sub_final <- dplyr::left_join(x_sub_filtered,dd,"id")
    x_sub_final <- x_sub_final[which(!is.na(x_sub_final$cluster)),]
    x_sub_final$cluster <- factor(as.character(x_sub_final$cluster))
    x_sub_final$Classification <- factor(x_sub_final$Classification,levels = unique(x_sub_final$Classification))
    x_sub_final$Solvent <- factor(x_sub_final$Solvent,levels = unique(x_sub_final$Solvent))
    
    uni_clas <- unique(x_sub_final$Classification)

    
    lm2_list <- list()
    cor_list <- list()
    cor_table_list <- list()
    
    for(j in 1:length(uni_clas)){
      # print(j)
      sub_uni_class <- subset(x_sub_final,Classification %in% uni_clas[[j]])
      sub_uni_class$D <- sub_uni_class$D/100000
      if(nrow(sub_uni_class)>5){
      lm2_list[[j]] <- lm(sub_uni_class$D ~ sub_uni_class$Mm)
      cor_list[[j]] <- cor.test(sub_uni_class$D, sub_uni_class$Mm)
      
        cor_table_list[[j]] <- data.frame(
        name = uni_clas[[j]],
        solvent = as.character(Solvent_list[[i]]),
        cor_st = cor_list[[j]]$statistic,
        cor_p_value = cor_list[[j]]$p.value,
        cor_est = cor_list[[j]]$estimate,
        lm_interc = lm2_list[[j]]$coefficients[[1]],
        lm_x = lm2_list[[j]]$coefficients[[2]],
        sd_lm_interc = summary(lm2_list[[j]])$coefficients[1,2],
        sd_lm_x = summary(lm2_list[[j]])$coefficients[2,2],
        p_lm_interc=summary(lm2_list[[j]])$coefficients[1,4],
        p_lm_x=summary(lm2_list[[j]])$coefficients[2,4],
        r.squared=summary(lm2_list[[j]])$r.squared,
        r.squared.adj=summary(lm2_list[[j]])$adj.r.squared
        )
      } else {
        cor_table_list[[j]] <- data.frame(
          name = uni_clas[[j]],
          solvent = as.character(Solvent_list[[i]]),
          cor_st = NA,
          cor_p_value = NA,
          cor_est = NA,
          lm_interc = NA,
          lm_x = NA,
          sd_lm_interc=NA,
          sd_lm_x = NA,
          p_lm_interc=NA,
          p_lm_x=NA,
          r.squared=NA,
          r.squared.adj=NA
          ) 
        
        lm2_list[[j]] <- NULL
        cor_list[[j]] <- NULL
      }
      
    }
    
    cor_table_list <- do.call(rbind,cor_table_list)
    
   
    write.csv(cor_table_list, paste0(out_dir,"/","fix_20190930","/","COR_SOLUTE_",
                                     as.character(Solvent_list[[i]]),"_",
                     Sys.Date(),".csv"),row.names = F,quote = F)
    
    
    # summary(lm2)
    # anova(lm2)
    
    x_sub_final$color_1 <- NA
    
    colors_filtered <- colors
    colors_filtered <- colors_filtered[which(colors_filtered$classification %in% uni_clas),]
    
    
    for(k in 1:nrow(colors_filtered)){
      x_sub_final$color_1[which(x_sub_final$Classification %in% colors_filtered$classification[[k]])] <- as.character(colors_filtered$color[[k]])
    };rm(k)
    
    x_sub_final$Classification <- factor(x_sub_final$Classification,levels=colors_filtered$classification)
    
    # #plot(x_sub_final$Mm,x_sub_final$D,col=x_sub_final$cluster)
      #  x_graph <- ggscatter(x_sub, x = "Mm", y = "D",
    x_graph <- ggscatter(x_sub_final, x = "Mm", y = "D",
                          #colour="color_1",
                         size = 20,
                         alpha = 0.5,
                         fill =  "Classification",
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
                         
                         xlim=c(min(x_sub_filtered$Mm,na.rm=T)-10,max(x_sub_filtered$Mm,na.rm=T)+10),
                         ylim=c(min(x_sub_filtered$D,na.rm=T),max(x_sub_filtered$D,na.rm=T)),
                         # ylab = expression(paste(Dx10^{-5},"/cm"^{2},"S"^{-1})),#cm^{2}s0^{-1}
                         ylab = expression(paste("Diffusion coefficient (cm"^{2},"S"^{-1},") x 10"^{-5})),#cm^{2}s0^{-1} 
                         xlab = expression(paste("Molecular weight "," ", "(g mol"^{-1},")")),
                         conf.int=T,add = "reg.line", repel = TRUE,
                         color = "Classification",
                         palette = as.character(colors_filtered$color)#"Dark2"
                           ) +
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
      theme(panel.background = element_rect(fill = "white"),
            text=element_text(size=70),
            axis.text.x  =element_text(size=70,angle = 90),
            axis.text.y  = element_text(size=70,colour="black"),
            legend.title=element_text(size=70,colour="black"),
           axis.line = element_line(color = "black",size = 1.2,linetype = "solid"))+
      scale_x_continuous(breaks = seq(min(x_sub_filtered$Mm,na.rm=T), max(x_sub_filtered$Mm,na.rm=T)+40, 20),labels=scaleFUN)+ 
      scale_y_continuous(breaks = seq(min(x_sub_filtered$D,na.rm=T), max(x_sub_filtered$D,na.rm=T), y_breaks),labels=scaleFUN_Y) #+
  #   theme_light()
    #+
    # xlim(c(min(x_sub$Mm),max(x_sub$Mm)))#+
    # ylim(c(min(x_sub$D),max(x_sub$D)))
    #stat_cor(aes(color = "black"), method = "pearson")#+
    #rremove("legend")+
    #grids(linetype = "solid",color = "gray94",size = 1)
    
    #x_graph
    ggsave(paste0(out_dir,"/","fix_20190930","/MANUAL/",as.character(Solvent_list[[i]]),"_",
                  Sys.Date(),".png"),x_graph,units="in",width=22,height=19,scale=2,dpi=200)
    
    #########################SUGGESTED
    #http://uc-r.github.io/kmeans_clustering
    x_graph2 <- ggscatter(x_sub_final, x = "Mm", y = "D",
                         
                         size = 20,
                         alpha = 0.5,
                         fill =  "cluster",
                         ggtheme = theme_light(),
                         title = Solvent_list[[i]],
                         rug = F,# Add marginal rug
                         #facet.by="Solvent",
                         short.panel.labs=F,
                         #dot.size = 800,
                         legend.title="Compounds by K-means",
                         cor.coef = F,
                         
                         label = "Solute",
                         font.label=c(48, "plain"),
                         
                         xlim=c(min(x_sub_final$Mm)-10,max(x_sub_final$Mm)+10),
                         ylim=c(min(x_sub_final$D),max(x_sub_final$D)),
                         # ylab = expression(paste(Dx10^{-5},"/cm"^{2},"S"^{-1})),#cm^{2}s0^{-1}
                         ylab = expression(paste("Diffusion coefficient (cm"^{2},"S"^{-1},") x 10"^{-5})),#cm^{2}s0^{-1} 
                         xlab = expression(paste("Molecular weight "," ", "(g mol"^{-1},")")),
                         conf.int=T,add = "reg.line", repel = TRUE,
                         color = "cluster", palette = "Dark2") +
      stat_cor(
        aes(color=cluster,label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), 
        method = "pearson",size=20,
        label.x.npc=0.7,
        label.y.npc=0.95,#,#,
        # label.x = 3
      )+# +
      #stat_regline_equation(aes(color=Classification))
      grids(axis = "xy",linetype = "dashed")
    
    x_graph2 <- x_graph2 +
      bgcolor("#FFFFFF")+
      theme(panel.background = element_rect(fill = "white"),
            text=element_text(size=70),
            axis.text.x  =element_text(size=50,angle = 90),
            axis.text.y  = element_text(size=70,colour="black"),
            legend.title=element_text(size=70,colour="black"),
            axis.line = element_line(color = "black",size = 1,linetype = "solid"))+
      scale_x_continuous(breaks = seq(min(x_sub_filtered$Mm,na.rm=T), max(x_sub_filtered$Mm,na.rm=T)+40, 20),labels=scaleFUN)+ 
      scale_y_continuous(breaks = seq(min(x_sub_filtered$D,na.rm=T), max(x_sub_filtered$D,na.rm=T), y_breaks),labels=scaleFUN_Y) #+
    ggsave(paste0(out_dir,"/","fix_20190930","/CLUSTER/",as.character(Solvent_list[[i]]),"_CLUSTER_",
                  Sys.Date(),".png"),x_graph2,units="in",width=22,height=19,scale=2,dpi=200)
    
    } else {
      cat(paste0("NO RECORDS FOR: ",Solvent_list[[i]]),"\n")
      x_sub_final <- NULL
      
    }
    return(x_sub_final)
  })
  #https://www.datanovia.com/en/lessons/k-means-clustering-in-r-algorith-and-practical-examples/#targetText=K%2Dmeans%20algorithm%20requires%20users,of%20clusters%20to%20be%20produced.
  ####################################
  sub_filtered <- do.call(rbind,sub_filtered)
  sub_filtered <- as.data.frame(sub_filtered[which(sub_filtered$Solvent!="Hydroquinone"),])
  ####################################
  dat_to_graph_original2 <- sub_filtered
 # dat_to_graph_original2$D <- dat_to_graph_original2$D/100000
  order <- data.frame(
    solvent= names(tapply(dat_to_graph_original2$D,dat_to_graph_original2$Solvent,mean)),
    value= tapply(dat_to_graph_original2$D,dat_to_graph_original2$Solvent,mean)
  )
  
  order <- order[complete.cases(order),]
  order <-  order[order(order$value,decreasing=T),]
  dat_to_graph_original2$Solvent <- factor(dat_to_graph_original2$Solvent,levels = order$solvent)
  
  
  
  
  #dat_to_graph_original2$Mm <- sub_filtered$Mm#/100000
  
  dat_to_graph_original3 <- dat_to_graph_original2[which(dat_to_graph_original2$Solvent!="Dimethyl sulfoxide"),]
 
  order <- data.frame(
    solvent= names(tapply(dat_to_graph_original3$D,dat_to_graph_original3$Solvent,mean)),
    value= tapply(dat_to_graph_original3$D,dat_to_graph_original3$Solvent,mean)
    )
  
  order <- order[complete.cases(order),]
  order <-  order[order(order$value,decreasing=T),]
  dat_to_graph_original3$Solvent <- factor(dat_to_graph_original3$Solvent,levels = order$solvent)
 
  
  
   x_graph_box <- ggboxplot(dat_to_graph_original2, x = "Solvent", y = "D", size = 2,
                           #fill =  "Classification",
                           rug = F,# Add marginal rug
                           #facet.by="Solvent",
                           short.panel.labs=F,
                           
                           # ylab = expression(paste(Dx10^{-5},"/cm"^{2},"S"^{-1})),#cm^{2}s0^{-1}
                           ylab = expression(paste("D (cm"^{2},"S"^{-1},") x 10"^{-5})),#cm^{2}s0^{-1} 
                           xlab = "",
                           #conf.int=T,add = "reg.line",
                           repel = TRUE,
                           #color = "black", 
                           palette = "Dark2") +
    stat_cor(aes(color = "black"), method = "spearman")+
    # rremove("legend")+
    grids(linetype = "solid",color = "gray94",size = 1)+
    theme(axis.text.x = element_text(angle = 90,size = 30),
          axis.text.y = element_text(size = 30),
          axis.title=element_text(size=25,face="bold")
          )
  
  # x_graph_box
  
  
  x_graph_box
  
  x_graph_box2 <- ggboxplot(dat_to_graph_original2[which(dat_to_graph_original2$Solvent!="Dimethyl sulfoxide"),], x = "Solvent", y = "D", size = 2,
                           #fill =  "Classification",
                           rug = F,# Add marginal rug
                           #facet.by="Solvent",
                           short.panel.labs=F,
                           
                           # ylab = expression(paste(Dx10^{-5},"/cm"^{2},"S"^{-1})),#cm^{2}s0^{-1}
                           ylab = expression(paste("D (cm"^{2},"S"^{-1},") x 10"^{-5})),#cm^{2}s0^{-1} 
                           xlab = "",
                           #conf.int=T,add = "reg.line",
                           repel = TRUE,
                           #color = "black", 
                           palette = "Dark2") +
    stat_cor(aes(color = "black"), method = "spearman")+
    # rremove("legend")+
    grids(linetype = "solid",color = "gray94",size = 1)+
    theme(axis.text.x = element_text(angle = 90,size = 30),
          axis.text.y = element_text(size = 30),
          axis.title=element_text(size=25,face="bold")
    )
  
  # x_graph_box
  
  
  x_graph_box2 
  
  # sub_filtered$D <- sub_filtered$D/100000
  
  
  
  
  
  # openxlsx::write.xlsx(sub_filtered,paste0("E:/JAVERIANA/ROTACION","/","sub_filtered_mm.xlsx"))
  
