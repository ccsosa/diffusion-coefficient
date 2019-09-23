# install.packages("webchem")
require(pacman)
#Calling up libraries
pacman::p_load(webchem,googlesheets,ggplot2,dplyr,BBmisc,rcurl,rvest,classyfireR)#Chemmine
library(gtools)

# if (!requireNamespace("BiocManager", quietly=TRUE))
#   install.packages("BiocManager")
# BiocManager::install("ChemmineR")

#Getting Google sheets list in a personal Gmail account
(my_sheets <- gs_ls())
my_sheets %>% glimpse()

##Calling up the google sheet table desired

Sheet <- gs_title("Diffussion_coefficient_2019_09_02")
gs_ws_ls(Sheet)
#Reading table
dat <- gs_read(ss = Sheet, ws = "data")

#Getting solute list from solute column
solute_list <- unique(dat$Solute)

#Extracting molecular weight via  Webchem
x_mW <- lapply(1:length(solute_list),function(i){
  
  
  # if(BBmisc::is.error(cir_query(as.character(solute_list[[i]]), representation = 'mw',
  #                               ))){
  #   x <- NULL
  # } else {
    x <- suppressAll(cir_query(as.character(solute_list[[i]]), representation = 'mw' #,
                   ) #resolver="name_pattern")#resolver = c("name_pattern","name_by_chemspider","name_by_cir","name_by_opsin","smiles"))
    )
                    # }
    cat(round((i/length(solute_list)*100),3)," DONE!","\n")
  return(x)
  
})


# x_mW[1:100]

x_mW2 <- x_mW



x_mW2M <- lapply(1:length(x_mW),function(i){
  
  x <- x_mW[[i]][[1]][1]
  
  return(x)
})


# Saving Molecular weight output to be look it up
#
x_mW2M <- do.call(rbind,x_mW2M)
output <- data.frame(name=solute_list,mW= x_mW2M)
write.table(output,"E:/JAVERIANA/ASIGNATURAS/2019_2/ROTACION_DROCHSS/20190911/mW.csv",
            quote=F,sep = "|",na = "",row.names = F)



#Getting InChIKey to get compounds classification
x_class <- list()
#x_class <- lapply(1:length(solute_list),function(i){
for(i in 1:length(solute_list)){  
cat(i,"\ proccessing","\n")
x2 <- webchem::cir_query(as.character(solute_list[[i]]),"stdinchikey") #,
x_class[[i]] <- x2[[1]][1]   
}

#Getting a table with compounds name and InChIKey without NA´s

x_class_tab <- data.frame(name=solute_list,InChIKey=NA)
x_class_tab$InChIKey <- x_class
x_class_tab$InChIKey <- sub("InChIKey=","",x_class_tab$InChIKey)
x_class_tab <- x_class_tab[which(x_class_tab$InChIKey!="NA"),]

###Splitting up compounds in chunks of n (5) compounds to get classification
n <- 5
nr <- length(x_class_tab$InChIKey)
server.species_splits <- split(x_class_tab$InChIKey, rep(1:ceiling(nr/n), each=n, length.out=nr))


##Getting compounds classification using lags of 2 and 20 second each at to avoid be banned!
#RUn each time you are banned by the server until get your classification done!s
families <- list()

##for(j in 1:length(server.species_splits)){ #1 DONE!
#for(j in 11:length(server.species_splits)){ #2 DONE!
#for(j in 21:length(server.species_splits)){ #3 DONE!
## for(j in 31:length(server.species_splits)){ #4 DONE!
##for(j in 41:length(server.species_splits)){ #5 DONE!
 for(j in 51:length(server.species_splits)){ #6 DONE!

            
#  lapply(1:length(server.species_splits),function(j){
  cat(paste(j,"chunks of ",length(server.species_splits),"starting!"),"\n")
  Sys.sleep(2);cat("2 seconds lag to start","\n")
  classification_list <- purrr::map(server.species_splits[[j]], get_classification)
  Sys.sleep(20);cat("20 seconds lag to continue","\n")
  families[[j]] <- classification_list
  #return(classification_list)
  cat(paste(j,"chunks of ",length(server.species_splits),"processed!"),"\n")
}  
#})

#First round of 12 chunks

# families1 <- families #DONE!; 
# families1 <- families1[which(!sapply(families1, is.null))]
# 
# families2 <- families #DONE!
# families2 <- families2[which(!sapply(families2, is.null))]

# families3 <- families #DONE!
# families3 <- families3[which(!sapply(families3, is.null))]

# families4 <- families #DONE!
# families4 <- families4[which(!sapply(families4, is.null))]

 # families5 <- families #DONE!
 # families5 <- families5[which(!sapply(families5, is.null))]
 
# families6 <- families #DONE!
# families6 <- families6[which(!sapply(families6, is.null))]

# rm(families)
 

families_final <- c(families1,families2,families3,families4,families5,families6) 


compound_list <- lapply(1:length(families_final),function(i){
  cat(paste("Chunk: ",i),"\n")
n <- length(families_final[[i]])

  x_list <- lapply(1:n,function(j){
  # cat(j,"\n")
    cat(paste("Chunk: ",i," | "," j:",j),"\n")
  x <- data.frame(kingdom=NA,superclass=NA,class=NA,subclass=NA,level_5=NA,level_6=NA)
  x2 <- t(as.data.frame(families_final[[i]][[j]]))
  if(nrow(x2)>0){
  colnames(x2) <- x2[1,]
  x2 <- x2[2,];x2 <- as.data.frame(t(x2))
  if(is.null(x2$kingdom)){x$kingdom <- NA} else {x$kingdom <- x2$kingdom}
  if(is.null(x2$superclass)){x$superclass <- NA} else {x$superclass <- x2$superclass}
  if(is.null(x2$class)){x$class <- NA} else {x$class <- x2$class}
  if(is.null(x2$`level 5`)){x$level_5 <- NA} else {x$level_5 <- x2$`level 5`}
  if(is.null(x2$`level 6`)){x$level_6 <- NA} else {x$level_6 <- x2$`level 6`}
  } else {
    x <- data.frame(kingdom=NA,superclass=NA,class=NA,subclass=NA,level_5=NA,level_6=NA)
  }
    return(x)
                                  })
  x_list <- do.call(rbind,x_list)
return(x_list)
  })

compound_list <- do.call(rbind,compound_list)
compound_list <- cbind(x_class_tab,compound_list)
write.table(compound_list,"E:/JAVERIANA/ASIGNATURAS/2019_2/ROTACION_DROCHSS/20190911/compound.csv",
            quote=F,sep = "|",na = "",row.names = F)
