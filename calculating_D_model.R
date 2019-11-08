pacman::p_load(googlesheets,dplyr,Metrics,reshape2,ggpubr,ggplot2,cluster,factoextra,BBmisc,readxl,openxlsx)#Chemmine

inDir <- "E:/JAVERIANA/ROTACION"
source("E:/JAVERIANA/ROTACION/diffusion-coefficient-master/ECUACIONES.R")

data_to <- readxl::read_excel(paste0(inDir,"/","sub_filtered_mm.xlsx"),sheet = "Sheet 1")
data_to <- as.data.frame(data_to)
data_to$D_CALC <- NA

#parameters <- readxl::read_excel(paste0(inDir,"/","lm_summary_FINAl.xlsx"),sheet = "Sheet 1")
parameters <- readxl::read_excel(paste0(inDir,"/","lm_summary_FINAl.xlsx"),sheet = "Sheet 1")

parameters <- as.data.frame(parameters)
parameters$rmse <- NA

x_sub_par <- list()

for(i in 1:nrow(parameters)){
# 
# #Getting Google sheets list in a personal Gmail account
# (my_sheets <- gs_ls())
# my_sheets %>% glimpse()
# Sheet <- gs_title("Diffussion_coefficient_2019_09_02")
# gs_ws_ls(Sheet)
# #Reading table
# dat <- gs_read(ss = Sheet, ws = "data")




#func_1_16(k,temp,Dens,n,visc,M)
  

sub1 <- data_to[which(data_to$Classification==parameters$`Compound species`[[i]] &
                data_to$Solvent==parameters$solvent[[i]]),]

for(j in 1:nrow(sub1)){
sub1$D_CALC[[j]] <- func_1_16(k=k,
                temp=parameters$`Temp (K)`[[i]],
                Dens=sub1$Molecular_density[[j]],
                n=parameters$n_or[[i]],
                visc=parameters$`viscocity (Centipoises)`[[i]],
                M=parameters$`Mm0 (g/mol)`[[i]],
                Mw=sub1$Mm[[j]])
};rm(j)

parameters$rmse[[i]] <- rmse(sub1$D,sub1$D_CALC)

x_sub_par[[i]] <- sub1
};rm(i)


x_sub_par <- do.call(rbind, x_sub_par)

#openxlsx::write.xlsx(solv_x_par,paste0("E:/JAVERIANA/ROTACION","/","solv_x_par.xlsx"))
openxlsx::write.xlsx(x_sub_par,paste0("E:/JAVERIANA/ROTACION","/","sub_filtered_D_CALC.xlsx"))
openxlsx::write.xlsx(parameters,paste0("E:/JAVERIANA/ROTACION","/","parameters_rmse.xlsx"))

