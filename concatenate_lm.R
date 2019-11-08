pacman::p_load(googlesheets,dplyr,reshape2,ggpubr,ggplot2,cluster,factoextra,BBmisc,readxl,openxlsx)#Chemmine
###########################################################
isEmpty <- function(x) {
  return(length(x)==0)
}

decimalplaces <- function(x) {
  if ((x %% 1) != 0) {
    nchar(strsplit(sub('0+$', '', as.character(x)), ".", fixed=TRUE)[[1]][[2]])
  } else {
    return(0)
  }
}
###########################################################
#Getting Google sheets list in a personal Gmail account
(my_sheets <- gs_ls())
my_sheets %>% glimpse()
Sheet <- gs_title("Diffussion_coefficient_2019_09_02")
gs_ws_ls(Sheet)
#Reading table
dat <- gs_read(ss = Sheet, ws = "data")
colnames(dat)[1] <- "id"
dat$id <- as.character(dat$id)
###########################################################
InDir <-  "E:/JAVERIANA/ROTACION/Graficos/fix_20190930"

source("E:/JAVERIANA/ROTACION/diffusion-coefficient-master/ECUACIONES.R")

fil_to <- list.files(InDir,pattern=".csv$")

files <- lapply(1:length(fil_to),function(i){
  x <- read.csv(paste0(InDir,"/",fil_to[[i]]))
  return(x)
})


files <- do.call(rbind,files)
files <- files[which(!is.na(files$cor_est)),]
files$Mm0 <- NA
files$Mm0 <- abs((files$lm_interc/files$lm_x)/4)

# files$sutherland <- NA
openxlsx::write.xlsx(files,paste0("E:/JAVERIANA/ROTACION","/","lm_summary.xlsx"))

###########################################################

sub_filtered <- openxlsx::read.xlsx(paste0("E:/JAVERIANA/ROTACION","/","sub_filtered_mm.xlsx"))

sub_filtered2 <- dplyr::right_join(sub_filtered,dat,"id")
sub_filtered2 <- sub_filtered2[which(sub_filtered2$id %in% sub_filtered$id),]

solvents <- unique(sub_filtered2$Solvent.x)

solv_x_par <- lapply(1:length(solvents),function(i){
  cat(paste("solvent: ",i),"\n")
  x_sub <- subset(sub_filtered2,Solvent.x %in% solvents[[i]])
  classes <- unique(x_sub$Classification)
  
  temp_median <- median(x_sub$`Temp (K)`,na.rm = T)
  density_median <- median(x_sub$Molecular_density.y,na.rm = T)
  visc_median <- median(x_sub$Viscocity,na.rm = T)
  # n_a(k=k,temp=temp_median,D=density_median,visc=visc_median,M,a)
  
  sub_file <- subset(files,solvent %in% solvents[[i]])
  
  x_par <- lapply(1:length(classes),function(j){
    cat(paste("class: ",j),"\n")
    Mm0 <- files[which(files$solvent== solvents[[i]] & files$name== classes[[j]]),]$Mm0
    if(isEmpty(Mm0)){Mm0 <- NA}
    a_m <- files[which(files$solvent== solvents[[i]] & files$name== classes[[j]]),]$lm_interc #inter
    if(isEmpty(a_m)){a_m <- NA}
    #a_m <- a_m*1E5
       # a_m <- sprintf("%.5f", a_m);a_m <-  as.numeric(a_m)
    m_m <- files[which(files$solvent== solvents[[i]] & files$name== classes[[j]]),]$lm_x #slope
    if(isEmpty(m_m)){m_m <- NA}
    #m_m <- m_m*1E8
    n1 <- n_a_alt(k=k,temp=temp_median,D=density_median,visc=visc_median,M=Mm0,a=a_m)
    if(isEmpty(n1)){n1 <- NA}
    n2 <- abs(n_m_alt(k=k,temp=temp_median,D=density_median,visc=visc_median,M=Mm0,m=m_m))
    if(isEmpty(n2)){n2 <- NA}
    if(isEmpty(n1) & isEmpty(n2) | is.na(n1) & is.na(n2)){
      n_F <- NA
      n_F2 <- n_F
      } else {
    n_F <- mean(c(n1,n2),na.rm = T);n_F2 <- n_F
    # format(n_F, scientific=F)
    n_F <- n_F*10^abs(as.numeric(strsplit(as.character(n_F),split = "e")[[1]][2]))
    }
    
    
    x_par <- data.frame(i=i,
                        j=j,
                        classes=classes[[j]],
                        a=a_m,
                        m=m_m,
                        temp=temp_median,
                        dens = density_median,
                        visc= visc_median,
                        solvent=solvents[[i]],
                        Mm0=Mm0,
                        n=n_F,
                        n_or=n_F2,
                        n_round=round(n_F,0),
                        n_ceiling=ceiling(n_F),
                        n1=n1,
                        n2=n2)
    return(x_par)
    })

  x_par <- do.call(rbind,x_par)
  return(x_par)
})


solv_x_par <- do.call(rbind,solv_x_par)
solv_x_par <- solv_x_par[complete.cases(solv_x_par),]

openxlsx::write.xlsx(solv_x_par,paste0("E:/JAVERIANA/ROTACION","/","solv_x_par.xlsx"))
