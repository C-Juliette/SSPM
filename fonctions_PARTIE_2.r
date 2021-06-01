source("fonctions_PARTIE_1.r", encoding = "UTF-8")

vecteur_directeur <- function(vecteur){
  Nb_premiers <- c(2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97)
  a <- abs(vecteur[1])
  b <- abs(vecteur[2])
  if (a != 0  | b != 0){
    if ((a != 0)  & (b != 0)){
      while( (a<1 | b<1) ){
        a <- a*10
        b <- b*10
      }
      for (i in Nb_premiers){
        while(a%%i == 0 & b%%i == 0 & (a != i | b != i)){
          a <- a/i
          b <- b/i
        }
      }
    }
    if (a == 0 | b == 0){
      a <- 0
      b <- 1
    }
    if (a == b){
      a <- 1 
      b <- 1
    }
  }

  vecteur_directeur <- c(a,b)
  return(vecteur_directeur)
}

covariance_theorique_Y <- function(sigma_Z, r, direction, echelle = 1){ 
  c <- 2*r+1 #cote de la fenetre
  i <- 0
  H1 <- c()
  H2 <- c()
  Distance_pixels <- c()
  vect_dir <- vecteur_directeur(direction)
  h <- sqrt(vect_dir%*%vect_dir)
  h1 <- vect_dir[1]
  h2 <- vect_dir[2]
  while (h1*i <= c & h2*i <= c){
    H1 <- append(H1,h1*i)
    H2 <- append(H2,h2*i)
    Distance_pixels <- append(Distance_pixels,h*i)
    i <- i+1
  }
  Covariance_theorique <- (sigma_Z**2)*(c**2 - (H1+H2)*c +H1*H2)/ c**4
  for (j in i:(i+5)){
    Distance_pixels <- append(Distance_pixels,h*j)
    Covariance_theorique <- append(Covariance_theorique,0)
  }
  Distance_km <- Distance_pixels*echelle
  result <- as.data.frame(cbind(Distance_pixels, Distance_km, Covariance_theorique))
  return(result)
}

construction_data_frame <- function(sigma_Z, les_rayons, les_directions, echelle = 1, vecteur_directeur = TRUE){
  datfr <- as.data.frame(setNames(replicate(4,numeric(0), simplify = F),c("Distance_pixels","Covariance_theorique","Direction","Rayon")))
  for (r in les_rayons){
    for (dir in les_directions){
      if (vecteur_directeur == TRUE){dir <- vecteur_directeur(dir)}
      df <- covariance_theorique_Y(sigma_Z, r, dir, echelle)
      df$Correlation_theorique <- df$Covariance_theorique/df$Covariance_theorique[1]
      df$Direction <- str_c(as.character(dir[1]), as.character(dir[2]), sep = " - ")
      #df$Nombre_pt <- (2*t+1)**2   A VOIR PLUS TARD
      df$Rayon_pixels <- r
      df$Rayon_km <- echelle*df$Rayon_pixels
      datfr <- rbind(datfr, df)
    }
  }  
  return(datfr)
}


graphique_covariance_theorique <- function(sigma_Z,les_rayons,les_directions,xlabs ="Distance en pixels", ylabs ="Covariance théorique",  x = "Distance_pixels", y = "Covariance_theorique", rayon = "Rayon_pixels",  echelle = 1,  max = "", maxy_sup ="", maxy_inf="", relier = TRUE,vecteur_directeur = TRUE){
  df <- construction_data_frame(sigma_Z, les_rayons, les_directions, echelle)
  ################COVARIANCE###################@
  if(y == "Covariance_theorique"){
  if(x== "Distance_pixels"){p <- ggplot(data = df, aes(x = Distance_pixels, y = Covariance_theorique))}
  if(x== "Distance_km"){p <- ggplot(data = df, aes(x = Distance_km, y = Covariance_theorique))}
  if (rayon == "Rayon_pixels"){p <- p + geom_point(aes(shape=Direction, col = factor(Rayon_pixels)))}
  if (rayon == "Rayon_km"){p <- p + geom_point(aes(shape=Direction, col = factor(Rayon_km)))}

  if(relier == TRUE){
    if (rayon == "Rayon_pixels"){p <- p + geom_line(aes(shape=Direction, col = factor(Rayon_pixels)))}
    if (rayon == "Rayon_km"){p <- p + geom_line(aes(shape=Direction, col = factor(Rayon_km)))}
  }
  }
  ###########CORRELATION###############
  if(y == "Correlation_theorique"){
    if(x== "Distance_pixels"){p <- ggplot(data = df, aes(x = Distance_pixels, y = Correlation_theorique))}
    if(x== "Distance_km"){p <- ggplot(data = df, aes(x = Distance_km, y = Correlation_theorique))}
    if (rayon == "Rayon_pixels"){p <- p + geom_point(aes(shape=Direction, col = factor(Rayon_pixels)))}
    if (rayon == "Rayon_km"){p <- p + geom_point(aes(shape=Direction, col = factor(Rayon_km)))}
    
    if(relier == TRUE){
      if (rayon == "Rayon_pixels"){p <- p + geom_line(aes(shape=Direction, col = factor(Rayon_pixels)))}
      if (rayon == "Rayon_km"){p <- p + geom_line(aes(shape=Direction, col = factor(Rayon_km)))}
    }
  }
  
  p <- p + 
    scale_color_viridis(discrete = TRUE, option = "B")+
    #scale_shape(solid = T)+
    labs(x = xlabs ,
         y = ylabs,
         shape = "Direction",
         color = rayon)+
    theme(
      legend.position = c(.95, .95),
      legend.justification = c("right", "top"),
      legend.box.just = "right",
      legend.margin = margin(6, 6, 6, 6)
    )
  if (max != ""){p <- p + scale_x_continuous(breaks= pretty_breaks()) + xlim(min=0, max = max)
    }
  if (maxy_sup != ""){p <- p +  ylim(min=maxy_inf, max = maxy_sup) }
  return(p)
}



################## FIN PARTIE THEORIQUE #####################

################## DEBUT PARTIE EMPIRIQUE #####################

covariance_empirique <- function(M, direction, echelle = 1){ 
  vect_dir <- vecteur_directeur(direction)
  Distance_pixels <- c()
  Covariance_empirique <- c()
  Correlation_empirique <- c()
  h <- sqrt(vect_dir%*%vect_dir)
  h1 <- vect_dir[1]
  h2 <- vect_dir[2]
  dim1 <- dim(M)[1]
  dim2 <- dim(M)[2]
  i <- 0
  while (h1*i <= (dim1-1) & h2*i <= (dim2-1)){
    ligne <- h1*i
    colonne <- h2*i
    Distance_pixels <- append(Distance_pixels,h*i)
    grille_arriere <- M[(ligne+1):dim1,(colonne+1):dim2]
    sigma_grille_arriere <- sqrt(variance(grille_arriere))
    grille_avant <- M[(1:(dim1-ligne)),1:(dim2-colonne)]
    sigma_grille_avant <- sqrt(variance(grille_avant))
    cov <- mean(grille_arriere*grille_avant)-mean(grille_arriere)*mean(grille_avant)
    cor <- cov/(sigma_grille_arriere*sigma_grille_avant)
    Covariance_empirique <- append(Covariance_empirique,cov)
    Correlation_empirique <- append(Correlation_empirique,cor)
    i <- i+1
  }
  Distance_km <- Distance_pixels*echelle
  result <- as.data.frame(cbind(Distance_pixels, Distance_km, Covariance_empirique, Correlation_empirique))
 # result$Correlation_empirique <- Covariance_empirique / 
  return(result)
}

construction_data_frame_emp <- function(Z, les_rayons, les_directions, echelle = 1, vecteur_directeur = TRUE){
  datfr <- as.data.frame(setNames(replicate(4,numeric(0), simplify = F),c("Distance_pixels","Distance_km","Covariance_empirique","Correlation_empirique")))
  for (r in les_rayons){
    for (dir in les_directions){
      if (vecteur_directeur == TRUE){dir <- vecteur_directeur(dir)}
      Y <- moy_gliss(Z,r)
      df <- covariance_empirique(Y, dir, echelle)
      df$Direction <- str_c(as.character(dir[1]), as.character(dir[2]), sep = " - ")
      df$Rayon_pixels <- r
      df$Rayon_km <- echelle*df$Rayon_pixels
      datfr <- rbind(datfr, df)
    }
  }  
  return(datfr)
}

graphique_covariance_empirique <- function(Z,les_rayons,les_directions,  ylabs = "", xlabs ="", x = "Distance_pixels", y = "Covariance_empirique", rayon = "Rayon_pixels",  echelle = 1,  max = "", maxy_sup = "", maxy_inf = "" , relier = TRUE,vecteur_directeur = TRUE){
  df <- construction_data_frame_emp(Z, les_rayons, les_directions, echelle)

  ################COVARIANCE###################@
  if(y == "Covariance_empirique"){
    if(x== "Distance_pixels"){p <- ggplot(data = df, aes(x = Distance_pixels, y = Covariance_empirique))}
    if(x== "Distance_km"){p <- ggplot(data = df, aes(x = Distance_km, y = Covariance_empirique))}
    if (rayon == "Rayon_pixels"){p <- p + geom_point(aes(shape=Direction, col = factor(Rayon_pixels)))}
    if (rayon == "Rayon_km"){p <- p + geom_point(aes(shape=Direction, col = factor(Rayon_km)))}
    
    if(relier == TRUE){
      if (rayon == "Rayon_pixels"){p <- p + geom_line(aes(shape=Direction, col = factor(Rayon_pixels)))}
      if (rayon == "Rayon_km"){p <- p + geom_line(aes(shape=Direction, col = factor(Rayon_km)))}
    }
  }
  ###########CORRELATION###############
  if(y == "Correlation_empirique"){
    if(x== "Distance_pixels"){p <- ggplot(data = df, aes(x = Distance_pixels, y = Correlation_empirique))}
    if(x== "Distance_km"){p <- ggplot(data = df, aes(x = Distance_km, y = Correlation_empirique))}
    if (rayon == "Rayon_pixels"){p <- p + geom_point(aes(shape=Direction, col = factor(Rayon_pixels)))}
    if (rayon == "Rayon_km"){p <- p + geom_point(aes(shape=Direction, col = factor(Rayon_km)))}
    
    if(relier == TRUE){
      if (rayon == "Rayon_pixels"){p <- p + geom_line(aes(shape=Direction, col = factor(Rayon_pixels)))}
      if (rayon == "Rayon_km"){p <- p + geom_line(aes(shape=Direction, col = factor(Rayon_km)))}
    }
  }
  
  p <- p + 
    scale_color_viridis(discrete = TRUE, option = "B")+
    #scale_shape(solid = T)+
    labs(x = xlabs ,
         y = ylabs,
         shape = "Direction",
         color = rayon)+
    theme(
      legend.position = c(.95, .95),
      legend.justification = c("right", "top"),
      legend.box.just = "right",
      legend.margin = margin(6, 6, 6, 6)
    )
  if (max != ""){p <- p +  xlim(min=0, max = max) }
  if (maxy_sup != ""){p <- p +  ylim(min=maxy_inf, max = maxy_sup) }
  return(p)
}

Comparaison_sans_moyenne_covariance <- function(Z, sigma_Z,les_rayons,les_directions, df,  x = "Distance_pixels", rayon = "Rayon_pixels",  echelle = 1,  max = "", maxy_sup = "",maxy_inf ="", relier = TRUE,vecteur_directeur = TRUE){
  p1 <- graphique_covariance_theorique(sigma_Z, les_rayons = les_rayons, x = x, y = "Covariance_theorique", rayon = rayon, les_directions = les_directions, echelle = echelle, max = max, maxy_sup = maxy_sup, maxy_inf=maxy_inf)+
    geom_point(df, aes(x,y))+
    geom_line(df, aes(x,y))
  p2 <- graphique_covariance_empirique(Z, les_rayons = les_rayons, les_directions = les_directions, x = x, y = "Covariance_empirique",  max = max, maxy_sup = maxy_sup, maxy_inf=maxy_inf)
p <- plot_grid(p1,p2, ncol=1, nrow=2)
return(p)
}

Comparaison_sans_moyenne_correlation <- function(Z, sigma_Z,les_rayons,les_directions,  x = "Distance_pixels", rayon = "Rayon_pixels",  echelle = 1,  max = "", maxy_sup = "",maxy_inf ="", relier = TRUE,vecteur_directeur = TRUE){
  p1 <- graphique_covariance_theorique(sigma_Z, les_rayons = les_rayons, x = x, y = "Correlation_theorique", rayon = rayon, les_directions = les_directions, echelle = echelle, max = max, maxy_sup = maxy_sup, maxy_inf=maxy_inf)
  p2 <- graphique_covariance_empirique(Z, les_rayons = les_rayons, les_directions = les_directions, x = x, y = "Correlation_empirique",  max = max, maxy_sup = maxy_sup, maxy_inf=maxy_inf)
  p <- plot_grid(p1,p2, ncol=1, nrow=2)
  return(p)
}


################### COMPARAISON à RAYONS et DIRECTION FIXES ##################
######################### PLOT VERT ROUGE COVARIANCE #####################

plrs_realisation <- function(r, direction, loi,  nbrealisation = 10, nbligne = 100, nbcolonne = 100, esperance_Z = 0, sigma_Z = 1, p = 0.5, lambda = 0.7,  echelle = 1){
  if(loi == "rbinom"){Z <- matrix(rnorm(nbligne*nbcolonne,1, p), nrow = nbligne)}
  if(loi == "rnorm"){Z <- matrix(rnorm(nbligne*nbcolonne,esperance_Z,sigma_Z), nrow = nbligne)}
  if(loi == "rpois"){Z <- matrix(rnorm(nbligne*nbcolonne,lambda), nrow = nbligne)}
  Y <- moy_gliss(Z,r)
  df <- covariance_empirique(Y, direction, echelle)[,1:3]
  setnames(df, old = 'Covariance_empirique', new = 'Covariance_empirique_1')
  for (i in 3:(nbrealisation+2)) {
    if(loi == "rbinom"){Z <- matrix(rnorm(nbligne*nbcolonne,1, p), nrow = nbligne)}
    if(loi == "rnorm"){Z <- matrix(rnorm(nbligne*nbcolonne,esperance_Z,sigma_Z), nrow = nbligne)}
    if(loi == "rpois"){Z <- matrix(rnorm(nbligne*nbcolonne,lambda), nrow = nbligne)}
    Y <- moy_gliss(Z,r)
    Covariance_empirique <- covariance_empirique(Y, direction, echelle)[,3]
    df <- cbind(df,Covariance_empirique)
    setnames(df, old = 'Covariance_empirique', new = str_c('Covariance_empirique_', (i-1), sep=""))
  }
  Moyenne <- rowMeans(df[,3:(nbrealisation+2)])
  df <- cbind(df, Moyenne)
  return(df)
}

Comparaison_covariances <- function(df_empirique, df_theorique, max = ""){
  p <- ggplot()+
    geom_point(data = df_theorique, aes(x = Distance_pixels, y = Covariance_theorique), col = 1) +
    geom_line(data = df_theorique, aes(x = Distance_pixels, y = Covariance_theorique), col = 1) 
   
  
  debut <- 1
  fin <- (dim(df_empirique)[2]-3)
  
  for(i in debut:fin){
    nom = str_c("Covariance_empirique_", as.character(i), sep = "")
    p <- p + geom_point(data = df_empirique, aes_string( x = "Distance_pixels", y = nom),  col = 3, size = 0.2)+
      geom_line(data = df_empirique, aes_string( x = "Distance_pixels", y = nom),  col = 3, size = 0.2)
  }
  
  
  p <- p+ geom_point(data = df_empirique, aes( x = Distance_pixels, y = Moyenne), col = 2)+
    geom_point(data = df_empirique, aes( x = Distance_pixels, y = Moyenne), col = 2)+ylab("Covariance")+
    
  
  p <- p+ xlab("Distance en pixels")
  
  
  
  if (max != ""){p <- p +  xlim(min=0, max = max) }
  return(p)
}


######################### PLOT VERT ROUGE CORRELATION #################



plrs_realisation_cor <- function(r, direction, loi, nbrealisation = 10, nbligne = 100, nbcolonne = 100, esperance_Z = 0, sigma_Z = 1, p = 0.5, lambda = 0.7,  echelle = 1){
  if(loi == "rbinom"){Z <- matrix(rnorm(nbligne*nbcolonne,1, p), nrow = nbligne)}
  if(loi == "rnorm"){Z <- matrix(rnorm(nbligne*nbcolonne,esperance_Z,sigma_Z), nrow = nbligne)}
  if(loi == "rpois"){Z <- matrix(rnorm(nbligne*nbcolonne,lambda), nrow = nbligne)}
  Y <- moy_gliss(Z,r)
  df <- covariance_empirique(Y, direction, echelle)[,c(1,2,4)]
  setnames(df, old = 'Correlation_empirique', new = 'Correlation_empirique_1')
  for (i in 3:(nbrealisation+2)) {
    if(loi == "rbinom"){Z <- matrix(rnorm(nbligne*nbcolonne,1, p), nrow = nbligne)}
    if(loi == "rnorm"){Z <- matrix(rnorm(nbligne*nbcolonne,esperance_Z,sigma_Z), nrow = nbligne)}
    if(loi == "rpois"){Z <- matrix(rnorm(nbligne*nbcolonne,lambda), nrow = nbligne)}
    Y <- moy_gliss(Z,r)
    Correlation_empirique <- covariance_empirique(Y, direction, echelle)[,4]
    df <- cbind(df,Correlation_empirique)
    setnames(df, old = 'Correlation_empirique', new = str_c('Correlation_empirique_', (i-1), sep=""))
  }
  Moyenne <- rowMeans(df[,3:(nbrealisation+2)])
  df <- cbind(df, Moyenne)
  return(df)
}

Comparaison_correlation <- function(df_empirique, df_theorique, max = ""){
  p <- ggplot()+
    geom_point(data = df_theorique, aes(x = Distance_pixels, y = Correlation_theorique), col = 1) +
    geom_line(data = df_theorique, aes(x = Distance_pixels, y = Correlation_theorique), col = 1) 
  
  
  debut <- 1
  fin <- (dim(df_empirique)[2]-3)
  
  for(i in debut:fin){
    nom = str_c("Correlation_empirique_", as.character(i), sep = "")
    p <- p + geom_point(data = df_empirique, aes_string( x = "Distance_pixels", y = nom),  col = 3, size = 0.2)+
      geom_line(data = df_empirique, aes_string( x = "Distance_pixels", y = nom),  col = 3, size = 0.2)
  }
  
  
  p <- p+ geom_point(data = df_empirique, aes( x = Distance_pixels, y = Moyenne), col = 2)+
    geom_point(data = df_empirique, aes( x = Distance_pixels, y = Moyenne), col = 2)+ylab("Corrélation")
  
  
  
  if (max != ""){p <- p +  xlim(min=0, max = max) }
  return(p)
}






library(stringr)
library(ggplot2)
library(viridis)
library(cowplot)
library(scales)


Z <- matrix(rnorm(102*102, 0, 1), nrow = 100)
sigma_Z <- 1
les_rayons <- 2
les_directions <- list(c(0,1))



x = 0:100
y1 = c(0.04, 0.04*20/25, 0.04*15/25, 0.04*10/25,0.04*5/25)
y2 = rep(0, 96 )
y <-c(y1,y2)
y
length(y)

df <- cbind(x,y)
df <- as.data.frame(df)
df
graphique_covariance_empirique(Z, les_rayons = les_rayons, les_directions = les_directions, y = "Covariance_empirique",  max = 100)
  


les_rayons <- 2
les_directions <- list(c(0,1))

df1 <- construction_data_frame_emp(Z, les_rayons, les_directions, echelle= 1)
y1 <- c(0.04, 0.04*20/25, 0.04*15/25, 0.04*10/25,0.04*5/25)
y2 <- rep(0, dim(df1)[1]-5)
y <-c(y1,y2)
y
dim(df1)[1]
df1$th <- y
df <- cbind(df1, y)  
df1
head(df1)

ggplot(df1)+
  geom_line(aes(x = Distance_pixels, y = th), col = "blue")+
  geom_point(aes(x = Distance_pixels, y = th), col = "blue")+
  geom_line(aes(x = Distance_pixels, y = Covariance_empirique), col = 2)+
  geom_point(aes(x = Distance_pixels, y = Covariance_empirique), col = 2)+
  xlab("Distance between points")+
  ylab("Covariance")



