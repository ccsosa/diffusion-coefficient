#require(ggplot2);require(reshape2);require(tidyr);require(ggpubr)

out_dir <- "E:/JAVERIANA/ROTACION/Graficos"
source("E:/JAVERIANA/ROTACION/diffusion-coefficient-master/ECUACIONES.R")
pacman::p_load(googlesheets,dplyr,reshape2,ggpubr,ggplot2,cluster,factoextra,BBmisc,Metrics,readxl,openxlsx)#Chemmine

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
                                12,13,9,15)]# 87 is solvent
colnames(dat_to_graph_original) <- c("id","Solute","Classification","Mm","D",
                                     "Solvent","Viscocity","Molecular_density","temp")
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


sub_filtered <- lapply(1:length(Solvent_list),function(i){
  cat(paste0("i: ",i," | plot for: ",Solvent_list[[i]]),"\n")
  x_sub <- subset(dat_to_graph_original,Solvent %in% Solvent_list[[i]])
  #cx_sub$D <- x_sub$D*100000
  
  x_sub_filtered  <- x_sub
  # x_sub_filtered <- x_sub[which(x_sub$D <x_max & x_sub$D >x_min),]
  x_sub_filtered <- x_sub_filtered[which(x_sub$Mm <400 & x_sub$Mm >80),]
  x_sub_final <- x_sub_filtered
  
  for(j in 1:nrow(x_sub_final)){
    func_1_18(k=k,
              temp=x_sub_final$temp[[j]],
              D=x_sub_final$D[[j]],
              n=,visc,M)   
  }
 
})