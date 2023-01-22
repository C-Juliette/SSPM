library(ggplot2)
library(cowplot)
library(docstring)
library(stringr)
library(viridis)
library(data.table)
library(scales)
library(tidyverse)
#test

moy_gliss <- function(Z, r){ 
  nblignesZ <- dim(Z)[1]
  nbcolonnesZ <- dim(Z)[2]
  nblignesY <- nblignesZ - 2*r
  nbcolonnesY <- nbcolonnesZ - 2*r
  
  Y <- matrix(data = rep(0, nblignesY*nbcolonnesY), nrow = nblignesY)
  
  for(i in (r+1):(nblignesZ-r)){
    for(j in (r+1):(nbcolonnesZ-r)){
      fenetre <- Z[(i-r):(i+r), (j-r):(j+r)]
      Y[i-r,j-r] <- mean(fenetre)
    }
  }
  return (Y)
} 

variance <- function(v){
  return(mean(v**2, na.rm =T)- mean(v, na.rm =T)**2)
}



affichage_matrice <- function(M, r= "", paletteinf = "", palettesup = "", titre = "", nom_axeX = "", nom_axeY = "", echelle = "Echelle"){
  titre <- str_c(titre, " - Var = " ,as.character(round(variance(M), 4)),  sep = "")
  M <- reshape2::melt(M) # trois colonnes : les i, les j, les valeurs
  p <- ggplot(data = M, aes(x = Var1 -0.5 , y = Var2-0.5) )
  if(r != ""){p <- ggplot(M, aes(x = Var1 -0.5+r , y = Var2 -0.5+r, z= value, fill=value))}
  p <- p + 
    geom_tile(aes(fill = value))
    if(paletteinf != "" & palettesup != ""){p <- p+scale_fill_viridis_c(option = "B", direction = -1, limits = c(paletteinf,palettesup))}
    else{p <- p+scale_fill_viridis_c(option = "B", direction = -1)}
    p <- p +
    scale_y_reverse() +
    labs(title = titre,
         x = nom_axeX,
         y = nom_axeY) +
    guides(fill = guide_colorbar(title = echelle))+
    scale_y_continuous(breaks= pretty_breaks())+
    scale_x_continuous(breaks= pretty_breaks())
  return(p)
}



affichage_niveau <- function(M, r = "", titre = "Lissage de Y", nom_axeX = "", nom_axeY = "", echelle = "Echelle"){
  M <- melt(M)
  M <- as.data.frame(M)
  p <- ggplot(M, aes(x = Var1, y = Var2 - 0.5, z= value - O.5, fill=value))
    if(r != ""){p <- ggplot(M, aes(x = Var1 -0.5+r , y = Var2 -0.5+r, z= value, fill=value))}
   p <- p +
    geom_raster(interpolate=TRUE)+
    scale_fill_viridis_c(option = "B", direction = -1)+
    labs(title = titre,
         x = nom_axeX,
         y = nom_axeY) +
    guides(fill = guide_colorbar(title = echelle))+
    scale_y_continuous(breaks= pretty_breaks())+
    scale_x_continuous(breaks= pretty_breaks())
  return(p)
}

les_boxplots <- function(Z, Y, titre = "Boxplots", nom_axeX = "", nom_axeY = ""){
  Z <- as.vector(Z)
  Y <- as.vector(Y)
  dataZ <- cbind(Z, rep("Z", length(Z)))
  dataY <- cbind(Y, rep("Y", length(Y)))
  data <- as.data.frame(rbind(dataZ, dataY))
  setnames(data, old = c('V2', "Z"), new = c("Grille", "Valeurs"))
  data$Valeurs <- as.numeric(data$Valeurs)
  data$Grille <- as.factor(data$Grille)
  data$Grille <- factor(data$Grille , levels=c("Z", "Y"))
  
  p <- ggplot(data, aes(x = Grille, y = Valeurs, fill = Grille))+  
    geom_boxplot()+
    stat_summary(fun.y="mean", color = "blue")+
    labs(title = titre,
         x = nom_axeX,
         y = nom_axeY) 
  p
}

les_histogrammes <- function(Z,Y, titre = 'Histogrammes', nom_axeX = "", nom_axeY = ""){
  Z <- as.vector(Z)
  Z <- as.data.frame(Z)
  Y <- as.vector(Y)
  Y <- as.data.frame(Y)
  
  p <- ggplot() + 
    geom_histogram(data = Y, aes(x = Y, y=..count../sum(..count..)), color="black", fill="#01BDC2" )+
    geom_histogram(data = Z, aes(x = Z, y=..count../sum(..count..)), color="black", fill="#F8766D", alpha=0.5)+
    
    labs(title = titre,
         x = nom_axeX,
         y = nom_axeY) 
    
  return(p)
}


affichage_general <- function(Z, r, titre = "",paletteinf = "", palettesup = ""){
  Y <- moy_gliss(Z,r)
  p1 <- affichage_matrice(Z, titre =  "Grille Z",paletteinf = paletteinf, palettesup = palettesup) + coord_cartesian(xlim = c(0, dim(Z)[1]), ylim = c(0, dim(Z)[2]))
  p2 <- affichage_matrice(Y,  r, titre = "Grille Y",paletteinf = paletteinf, palettesup = palettesup)+ coord_cartesian(xlim = c(0, dim(Z)[1]), ylim = c(0, dim(Z)[2]))
  #p3 <- affichage_niveau(Y,r) + coord_cartesian(xlim = c(0, dim(Z)[1]), ylim = c(0, dim(Z)[2]))
  p3 <- les_histogrammes(Z,Y)
  p4 <- les_boxplots(Z,Y)
  
  title <- ggdraw() + draw_label(titre, fontface='bold')
  
  p <- plot_grid(p1,p2,p3,p4, ncol=2, nrow=2) #, labels=c(str_c("Grille Z - Var(Z) = ", as.character(round(variance(Z), 4)), sep = ""),str_c("Grille Y - Var(Y) = ", as.character(round(variance(Y), 4)), sep = ""),"Lissage de Y", "Boxplots") )+
    scale_fill_viridis_c(option = "B", direction = -1)
  #plot_grid(titre, p, ncol=1, rel_heights=c(0.1, 1))
  p
}


affichage_avec_r_qui_augmente<- function(Z){
   p1 <- affichage_matrice(Z, titre = "Grille Z")
   p2 <- affichage_matrice(moy_gliss(Z,r =1), r = 1, titre = str_c("r=1 :  ", as.character((2*1+1)**2), " px/fenêtre " , sep = "")) + coord_cartesian(xlim = c(0, dim(Z)[1]), ylim = c(0, dim(Z)[2]))
   p3 <- affichage_matrice(moy_gliss(Z,r =2), r = 2, titre = str_c("r = 2 :  ", as.character((2*2+1)**2), " px/fenêtre " , sep = "")) + coord_cartesian(xlim = c(0, dim(Z)[1]), ylim = c(0, dim(Z)[2]))
   p4 <- affichage_matrice(moy_gliss(Z,r =3), r = 3, titre = str_c("r = 3 :  ", as.character((2*3+1)**2), " px/fenêtre " , sep = "")) + coord_cartesian(xlim = c(0, dim(Z)[1]), ylim = c(0, dim(Z)[2]))
   p5 <- affichage_matrice(moy_gliss(Z,r =5), r = 5, titre = str_c("r = 5 :  ", as.character((2*5+1)**2), " px/fenêtre " , sep = "")) + coord_cartesian(xlim = c(0, dim(Z)[1]), ylim = c(0, dim(Z)[2]))
   p6 <- affichage_matrice(moy_gliss(Z,r =8), r = 8, titre = str_c("r = 8 :  ", as.character((2*8+1)**2), " px/fenêtre " , sep = "")) + coord_cartesian(xlim = c(0, dim(Z)[1]), ylim = c(0, dim(Z)[2]))
   p7 <- affichage_matrice(moy_gliss(Z,r =10), r = 10, titre = str_c("r = 10 :  ", as.character((2*10+1)**2), " px/fenêtre " , sep = "")) + coord_cartesian(xlim = c(0, dim(Z)[1]), ylim = c(0, dim(Z)[2]))
   p8 <- affichage_matrice(moy_gliss(Z,r = 15), r = 15, titre = str_c("r = 15 :  ", as.character((2*15+1)**2), " px/fenêtre " , sep = "")) + coord_cartesian(xlim = c(0, dim(Z)[1]), ylim = c(0, dim(Z)[2]))
   p9 <- affichage_matrice(moy_gliss(Z,r =20), r = 20, titre = str_c("r = 20 :  ", as.character((2*20+1)**2), " px/fenêtre " , sep = "")) + coord_cartesian(xlim = c(0, dim(Z)[1]), ylim = c(0, dim(Z)[2]))
   
 p <- plot_grid(p1, p2, p3, p4, p5, p6, p7, p8, p9, ncol=3, nrow = 3)
 return(p)
}


histo_data_Y_correles <- function(Z){
  dframe <-c() 
  for(r in 1:25){
    legende <- c()
    Rayon <- c()
    Nb_points_fenetre <- c()
    Nb_points_grille <- c()
    for (i in 1:((100-2*r)**2)){
      Rayon[i] <- r
      Nb_points_fenetre[i] <- (2*r+1)**2
      Nb_points_grille[i] <- (100-2*r)**2
      #legende[i] <- str_c("fafaz", as.character(i), sep = "")
    }
    Y <- as.vector(moy_gliss(Z, r))
    df <- cbind.data.frame(Y, Rayon, Nb_points_fenetre, Nb_points_grille)
    dframe <- rbind.data.frame(dframe, df)
    #Legende <- rbind(Legende, legende)
  }
  dframe$Legende <- str_c("F : ", (2*dframe$Rayon+1)^2, " px = (2*", dframe$Rayon,"+1)² = " , (2*dframe$Rayon+1), "²", sep ="")
  #print(length(Legende))
  return(dframe)
}

histo_data_Y_independants <- function(loi, p=0.5, esperance_Z = 0, sigma_Z = 1, lambda = 0.7){
  W <- c()
  Rayon <- c()
  dframe <-c() 
  for(r in 1:25){
    for (i in 1:10000){
      if(loi == "rbinom"){W[i] <- mean(rbinom((2*r+1)**2,1, p))}
      if(loi == "rnorm"){W[i] <- mean(rnorm((2*r+1)**2, esperance_Z, sigma_Z))}
      if(loi == "rpois"){W[i] <- mean(rpois((2*r+1)**2, lambda))}
      Rayon[i] <- r
    }
    df <- cbind(W, Rayon)
    dframe <- rbind(dframe, df)
    dframe <- as.data.frame(dframe)
  }
  return(dframe)
}






