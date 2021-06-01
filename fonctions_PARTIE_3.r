source("fonctions_PARTIE_1.r", encoding = "UTF-8")
source("fonctions_PARTIE_2.r", encoding = "UTF-8")

##### Partie 3 #####

maillage_large <- function(Y, facteur_dilatation, bord = 0){
  Mlarge <- matrix(rep(NA, dim(Y)[1]*dim(Y)[2]), nrow = dim(Y)[1])
  abscisse <- seq(bord+1,nrow(Y)-bord,facteur_dilatation)
  ordonnee <- seq(bord+1,ncol(Y)-bord,facteur_dilatation)
  Mlarge[abscisse, ordonnee] <- Y[abscisse, ordonnee]
  return(Mlarge)
}


maillage_large_sans_NA <- function(Y, facteur_dilatation, bord = 0){
  return(Y[seq(bord+1,nrow(M)-bord,facteur_dilatation),seq(bord+1,ncol(M)-bord,facteur_dilatation)])
}


maillage_large_par_moyenne <- function(Y, rprime, grand_pixel = TRUE){
  facteur_dilatation <- 2*rprime+1
  
  Ytemp <- moy_gliss(Y, rprime)
  selection_ligne <- seq(1, dim(Ytemp)[1], facteur_dilatation)
  selection_colonne <- seq(1, dim(Ytemp)[2], facteur_dilatation)
  Y2_sansNA <- Ytemp[selection_ligne,selection_colonne]
  
  if(grand_pixel == TRUE){
    Y2 <- matrix(rep(0, dim(Y)[1]*dim(Y)[2]), nrow = dim(Y)[1])
    for(i in 1:dim(Y2_sansNA)[1]){
      for(j in 1:dim(Y2_sansNA)[2]){
        fenetre_ligne <- (facteur_dilatation*i - (facteur_dilatation -1)):(facteur_dilatation*i)
        fenetre_colonne  <- (facteur_dilatation*j - (facteur_dilatation -1)):(facteur_dilatation*j)
        Y2[fenetre_ligne,fenetre_colonne] <- Y2_sansNA[i,j]
      }
    }
  }
  if(grand_pixel == FALSE){Y2 <- Y2_sansNA}
  return(Y2)
}


s2intra <- function(Y, rprime){
  facteur_dilatation <- 2*rprime+1
  V <- c()
  bloc_ligne <- dim(Y)[1] %/% facteur_dilatation
  bloc_colonne <- dim(Y)[2] %/% facteur_dilatation 
  for( i in 1:bloc_ligne){
    for(j in 1:bloc_colonne){
      fenetre_ligne <- (facteur_dilatation*i - (facteur_dilatation -1)):(facteur_dilatation*i)
      fenetre_colonne  <- (facteur_dilatation*j - (facteur_dilatation -1)):(facteur_dilatation*j)
      fenetre <- Y[fenetre_ligne,fenetre_colonne]
      V <- append(V,variance(as.vector(fenetre)))
    }
  }
  return(mean(V))
}


s2inter <- function(Y,rprime){
  Y2 <- maillage_large_par_moyenne(Y, rprime = 2, grand_pixel = FALSE)
  sinter <- variance(Y2)
  return(sinter)
}


