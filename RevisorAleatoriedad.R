#########################################################
#                                                     #
#     REVISIONES de los generadores "aleatorios"      #
#               de R, C y un cuantico (H)             #
#                                                     #
#######################################################
 
# David Cortés Servín
# Guillermo Arriaga García
# Edgar Josué Pedraza Cervantes
 
# Proyecto de Probabilidad, 2013.
 
########################################################
##########################################   Descripción
 
 
# Pruebas de Uniformidad:   X2, KS, AD, Monte Carlo
# Pruebas de Independencia: Corridas, Serial
 
# Se prueban 100 muestras
#    se ven los p-valores o estadísticos
# Los aleatorios cuanticos se consiguen en
#   https://www.fourmilab.ch/hotbits/secure_generate.html
 
 
  
 
  
########################################################
#################################################   Test
 
  
  
####################### ### ## # Prueba Anderson Darling
 
Prueba.AD.unif01 <- function(muestraunif) {
 
  # Requiere a la muestra en [0,1]
 
  muestra <- sort(muestraunif)
 
  # Arreglo de ceros y unos
  ajuste1 <- 0.01*muestra[which(muestra>0)[1]]
  ajuste2 <- 0.01*muestra[which(muestra>ajuste1)[1]]
  ajuste  <- (ajuste1+ajuste2+0)/3  #Promedio (subjetivo)
 
 
 
  muestra[which(muestra==0)]=ajuste
  muestra[which(muestra==1)]=1-ajuste
 
  N <- length(muestra)
  S <- 1:N
 
 
  # Estadistico A*A=-N-S
  for(i in 1:N){
    S[i]=((2*i-1)/N)*(log(muestra[i])+log(1-muestra[N+1-i])) 
  }
 
  s <- sum(S)
  A <- -N-s
 
  # Valores críticos del test Anderson Darling
  # con todos los parametros conocidos y N > 4
  #
  # valor  1-a (significancia)
  # 1.933  0.90
  # 2.492  0.95   (común)
  # 3.070  0.975
  # 3.857  0.99
 
  # Región crítica: Rechazar H0 si A*A > 2.492
  # H0 es que los datos sigan la distribucion especificada
 
  #if(A > 2.492)
  #  print(paste(
  #   "Se rechaza H0. A*A=",round(A,4),"> 2.495 (95%)"))
  #else
  #  print(paste(
  #   "No se rechaza H0. A*A=",round(A,4),"< 2.495 (95%)"))
 
   return (A)
}
 
 
 
#################################### ### ## #  Prueba X2
 
Prueba.X2.unif01 <- function(m1,
                    clases=as.integer(sqrt(length(m1)))){
 
  # Requiere a la muestra en [0,1]
 
  m2 <- sort(m1)
  n2 <- length(m2)
 
  obs <- 1:clases
 
 
 
    # Por esto raiz de n clases: 
  for(i in 1:clases){
     obs[i]=length(which(m2<(i/clases)))
  }
  obs <- c(obs[1],diff(obs))
 
 
  observed <- obs
  expected <- rep( n2/clases , clases )
 
  TE <- cbind(observed,expected)
 
  AA <- sum( (observed-expected)^2/expected )
  aT <- 1-pchisq(AA,df=(clases-1))
 
  return( aT ) 
}
 
 
########################## ### ## # Método de Montecarlo
 
 
   g<-function(x){
        return(x^2)
   }
   no_negativo<-function(dominio){
      return( min(g(dominio)) )
   }
   f<-function(x){
      return( g(x)-no_negativo((0.001)*(1:1000)) )
   }
 
Montecarlo <- function(muestra, particion=2*floor((sqrt(length(muestra))/2)), impresion=F){
  
   k<-floor(length(muestra)/2)
   P<-cbind(muestra[1:k],muestra[(k+1):(2*k)])

   # Estimacion con la funcion f  
   aux<-0
 
   for(i in 1:k){
      if(f(P[i,1])>=P[i,2]){ 
         aux<-aux+1
      }
   }
 
   aux1<-aux
 
   Q<-matrix(ncol=2,nrow=aux)
 
   aux<-0
   for(i in 1:k){
      if(f(P[i,1])>=P[i,2]){
         aux<-aux+1
         Q[aux,1]<-P[i,1]
         Q[aux,2]<-P[i,2]
      }
   }
 
   if(impresion){
      x11()
      par(mfrow=c(1,2))
      plot(Q[,1],Q[,2],col="red",type="p",xlim=c(0,1),ylim=c(0,1),
           main="Montecarlo en una función",
           xlab="",ylab=""
      )
  
      par(new=TRUE)
      X<-0.001*(1:1000)
      plot(X,f(X),type="l",lwd=2,col="blue")
   }


# Estimacion con el tablero de ajedrez

   n<-particion

   A <- matrix(0,n,n)
   eje <- seq(0,1,length=n )[-1]
 
 
   R<-n*P
   aux<-0
   for(i in 1:k){
      a<-floor(R[i,1])
      b<-floor(R[i,2])
      if( ((a+b)%%2 ) ==0){
         aux<-aux+1
##         A[a+1,b+1]<-A[a+1,b+1]+1
      }
    }
 
   S<-matrix(ncol=2,nrow=aux)
   R<-n*P
   aux<-0
   for(i in 1:k){
      a<-floor(R[i,1])
      b<-floor(R[i,2])
      if(((a+b)%%2)==0){
         aux<-aux+1
         S[aux,1]<-P[i,1]
         S[aux,2]<-P[i,2]
      }
   }
  
   if(impresion){
      plot(S[,1],S[,2],col="red",type="p",xlim=c(0,1),ylim=c(0,1),main="Montecarlo")
   }


   # Conteos en cada casilla negra ajustados por Poisson con ks.test
   #    Entre mejor sean explicados por Poisson, mas uniforme y aleatoria es la muestra
   pos_par   <- (1:(n/2))*2
   pos_impar <- pos_par-1
 
#   Cuad_Negros <- c()
#   for(i in pos_par ){
#      Cuad_Negros=c(Cuad_Negros,A[i,pos_par])
#   }
#   for(i in pos_impar ){
#      Cuad_Negros=c(Cuad_Negros,A[i,pos_impar])
#   }
# 
#   cn <- length(Cuad_Negros)
#   p.val <- ks.test(Cuad_Negros,qpois(1:cn/(cn+1), mean(Cuad_Negros)) )$p.value
# 
   return(c(aux1/k,aux/k,0.95))
}
 
 
validacion <- function(){
  
   # Para identificar el tamaño de partición
   #   por eje como raiz del tamaño de muestra
 
   # Los conteos por regiones negras tipo tablero
   #   de ajedrez es poisson si se usa raiz de n
 
   n <- 10
   k<-100
  
   P<-cbind(runif(k,0,1),runif(k,0,1))
   A <- matrix(0,n,n)
   eje <- seq(0,1,length=n )[-1]
 
   R<-n*P
  
   aux<-0
   for(i in 1:length(P[,1])){
      a<-floor(R[i,1])
      b<-floor(R[i,2])
 
      if(((a+b)%%2)==0){
         A[a+1,b+1]<-A[a+1,b+1]+1
      }
   }
  
   pos_par   <- (1:(n/2))*2
   pos_impar <- pos_par-1
 
   Cuad_Negros <- c()
   print( length( Cuad_Negros ) )
   for(i in pos_par ){
      Cuad_Negros=c(Cuad_Negros,A[i,pos_par])
      print( length( Cuad_Negros ) )
   }
   for(i in pos_impar ){
      Cuad_Negros=c(Cuad_Negros,A[i,pos_impar])
      print( length( Cuad_Negros ) )
   }
 
 
   plot( sort(Cuad_Negros) , lwd=3 , col=3)
   par(new=T)
 
   n2 <- length(Cuad_Negros)
 
   plot( sort(rpois(k,mean(Cuad_Negros))) , lwd=3 , col=4)
 
   plot( qpois(1:n2/(n2+1),mean(Cuad_Negros)) , sort(Cuad_Negros) , main=paste("Validacion QQ Poisson L=",round(mean(Cuad_Negros),3))  )
   for(i in 1:50){
      points( qpois(1:n2/(n2+1),mean(Cuad_Negros)), sort(rpois(n2,mean(Cuad_Negros))) , col=6)
   }
   lines( qpois(1:n2/(n2+1),mean(Cuad_Negros)) , sort(Cuad_Negros) )
}  





# PRUEBAS DE INDEPENDENCIA:
 

###############################  ## ### ## # Serial test

Prueba.Serial<-function(Muestra, k=floor(sqrt(length(Muestra)))){
   n <- length(Muestra)
   Mat <- matrix(0,k,k)
 
   # El [0,1] se parte en k=particion casillas
   # Se identifica cada dato en qué casilla cae

   posicion<-floor((Muestra*k))
   for(i in 2:n){
      Mat[posicion[i-1],posicion[i]]<-Mat[posicion[i-1],posicion[i]]+1
   }
   Sumas<-c(1:k)
   for(i in 1:k)
   Sumas[i]<-sum(Mat[i,])
   A<-((k^2)/n)*sum((Mat-(n/(k^2)))^2)-(k/n)*sum((Sumas-n/k)^2)

   p_val <- 1-pchisq(A,df=(k)*(k-1))
 
   return(p_val)
}


#################################  ## ### ## # Runs test
Prueba.Rachas<-function(Muestra){
   Rachas <- as.vector(sign(diff(Muestra)))
   Rachas[Rachas==0] <- 1
   n <- length(Muestra)-1

   for(i in 2:(n-1)){
      if(is.numeric(Rachas[i]) & is.numeric(Rachas[i-1]) ){
         if(Rachas[i]*Rachas[i-1]>0){
            Rachas[i]=Rachas[i-1]+Rachas[i]
            Rachas[i-1]=0
         }
      }
   }

   Rachas=sort(abs(Rachas[Rachas!=0]))
   Ra <-as.vector(as.numeric(levels(as.factor(Rachas))))
   Ra_obs <- Ra
   Ra_esp <- Ra
  
  
   for(k in Ra){
      Ra_obs[which(Ra==k)]=length(which(Rachas==k))
      Ra_esp[which(Ra==k)]=2*( n*(k*k+3*k+1)-(k^3+3*k*k-k-4) )/factorial(k+3)
   }
  
   A<-1:n
   for(i in 1:length(Ra)){
      A=A[-which(A==Ra[i])]
   }

   Estadistico <- sum(((Ra_esp-Ra_obs)^2)/Ra_esp)

   if(!is.na(which(A==n)[1])){
      A=A[A!=n]
      if(n<100){
         Estadistico <- Estadistico+2/factorial(n+1)
      }
   }

   for(k in A){
      Estadistico <- Estadistico+2*( n*(k*k+3*k+1)-(k^3+3*k*k-k-4) )/factorial(k+3)
   } # Valen cero estas longitudes "observadas como nulas", por eso la simplificación.

   p_val <- 1-pchisq(Estadistico,df=(k-1))

   return(p_val)
}

 
 
 
########################################################
###################################   Función de Pruebas
 
   # Recibe una muestra de 100,000 datos en (0,1)
  
 
Prueba <- function(muestra,tamano=c(100,1000),generador="",impresion=T){
 
   resAD  <- 1:100
   resKS  <- 1:100
   resX2  <- 1:100
   resMC1 <- 1:100
   resMC2 <- 1:100
   resMC3 <- 1:100
   resRT  <- 1:100
   resST  <- 1:100
#   resTC  <- 1:100
 
   RESULTADOS <- data.frame(
     cbind(Tamano=tamano,
           AD.90a=0, AD.95a=0, AD.99a=0,
           AD.90b=0, AD.95b=0, AD.99b=0,
           KS.90a=0, KS.95a=0, KS.99a=0,
           KS.90b=0, KS.95b=0, KS.99b=0,
           X2.90a=0, X2.95a=0, X2.99a=0,
           X2.90b=0, X2.95b=0, X2.99b=0,
           MC.Am=0,  MC.Av=0,
           MC.Bm=0,  MC.Bv=0,
           MC.Cm=0,  MC.Cv=0,
           RT.90a=0, RT.95a=0, RT.99a=0,
           RT.90b=0, RT.95b=0, RT.99b=0,
           ST.90a=0, ST.95a=0, ST.99a=0,
           ST.90b=0, ST.95b=0, ST.99b=0 #,
#           TC.90a=0, TC.95a=0, TC.99a=0,
#           TC.90b=0, TC.95b=0, TC.99b=0

     )
   )
   RESULTADOS_INDEP <- data.frame(
     cbind(Tamano=c(1000,10000,100000),
           RT.90a=0, RT.95a=0, RT.99a=0,
           RT.90b=0, RT.95b=0, RT.99b=0,
           ST.90a=0, ST.95a=0, ST.99a=0,
           ST.90b=0, ST.95b=0, ST.99b=0 #,
#           TC.90a=0, TC.95a=0, TC.99a=0,
#           TC.90b=0, TC.95b=0, TC.99b=0

     )
   )
  
 
 
   for(n in tamano ) {
 
      ## APLICACION DE LOS TEST
      for(i in 1:100){
         mm <- muestra[((i-1)*n+1):(n*i)]
         # Monte Carlo esta para uniformidad
         MCarlo <- as.vector(Montecarlo(mm) ) 
           resMC1[i]=MCarlo[1]
           resMC2[i]=MCarlo[2]
           resMC3[i]=MCarlo[3]
           resRT[i]=Prueba.Rachas(mm)
           resST[i]=Prueba.Serial(mm)
#           resTC[i]=Tablas.Contingencia(mm)
         mm <- sort(mm)
         resAD[i]=Prueba.AD.unif01(mm)
         resKS[i]=ks.test( mm, (1:n)/(n+1) )$p.value
         resX2[i]=Prueba.X2.unif01(mm)

      }
    
 
     resAD=sort(resAD)
     resKS=sort(resKS)
     resX2=sort(resX2)
     resMC1=sort(resMC1)
     resMC2=sort(resMC2)
     resMC3=sort(resMC3)
     resRT=sort(resRT)
     resST=sort(resST)
#     resTC=sort(resTC)
    
 
     
  # Cantidad de muestras no mayores a los niveles
  aAD <- (which(resAD>1.933)[1]-1)/100
  bAD <- (which(resAD>2.492)[1]-1)/100
  cAD <- (which(resAD>3.857)[1]-1)/100
 
  # Valores de rechazo estimados en los niveles
  aaAD <- resAD[90]
  bbAD <- resAD[95]
  ccAD <- resAD[99]
 
 
  aKS <- (100-(which(resKS>0.10)[1]-1))/100
  bKS <- (100-(which(resKS>0.05)[1]-1))/100
  cKS <- (100-(which(resKS>0.01)[1]-1))/100
 
  aaKS <- resKS[11] # Nivel de rechazo estimado al 90%
  bbKS <- resKS[6]  # ... al 5%
  ccKS <- resKS[2]  # ... al 1%
 
  aX2 <- (100-(which(resX2>0.10)[1]-1))/100
  bX2 <- (100-(which(resX2>0.05)[1]-1))/100
  cX2 <- (100-(which(resX2>0.01)[1]-1))/100
     
  aaX2 <- resX2[11] # Nivel de rechazo estimado al 90%
  bbX2 <- resX2[6]  # ... al 5%
  ccX2 <- resX2[2]  # ... al 1%
 
       
  vec <- correccion_NA(c(aAD,bAD,cAD))
  aAD <- vec[1]
  bAD <- vec[2] 
  cAD <- vec[3]
 
  vec <- correccion_NA(c(aKS,bKS,cKS))
  aKS <- vec[1]
  bKS <- vec[2] 
  cKS <- vec[3]
     
  vec <- correccion_NA(c(aX2,bX2,cX2))
  aX2 <- vec[1]
  bX2 <- vec[2] 
  cX2 <- vec[3]
     
 
 
  # Almacenamiento de resultados
  i <- which(tamano==n)[1]
  RESULTADOS$AD.90a[i]=aAD
  RESULTADOS$AD.95a[i]=bAD
  RESULTADOS$AD.99a[i]=cAD
  RESULTADOS$AD.90b[i]=aaAD
  RESULTADOS$AD.95b[i]=bbAD
  RESULTADOS$AD.99b[i]=ccAD
  RESULTADOS$KS.90a[i]=aKS
  RESULTADOS$KS.95a[i]=bKS
  RESULTADOS$KS.99a[i]=cKS
  RESULTADOS$KS.90b[i]=aaKS
  RESULTADOS$KS.95b[i]=bbKS
  RESULTADOS$KS.99b[i]=ccKS
  RESULTADOS$X2.90a[i]=aX2
  RESULTADOS$X2.95a[i]=bX2
  RESULTADOS$X2.99a[i]=cX2
  RESULTADOS$X2.90b[i]=aaX2
  RESULTADOS$X2.95b[i]=bbX2
  RESULTADOS$X2.99b[i]=ccX2
  RESULTADOS$MC.Am[i]=mean(resMC1)
  RESULTADOS$MC.Av[i]=var(resMC1)
  RESULTADOS$MC.Bv[i]=var(resMC2)
  RESULTADOS$MC.Bm[i]=mean(resMC2)
  RESULTADOS$MC.Cm[i]=mean(resMC3)
  RESULTADOS$MC.Cv[i]=var(resMC3)

  aRT <- (100-(which(resRT>0.10)[1]-1))/100
  bRT <- (100-(which(resRT>0.05)[1]-1))/100
  cRT <- (100-(which(resRT>0.01)[1]-1))/100
 
  aST <- (100-(which(resST>0.10)[1]-1))/100
  bST <- (100-(which(resST>0.05)[1]-1))/100
  cST <- (100-(which(resST>0.01)[1]-1))/100
 
  aaST <- resST[11] # Nivel de rechazo estimado al 90%
  bbST <- resST[6]  # ... al 5%
  ccST <- resST[2]  # ... al 1%

  aaRT <- resRT[11] # Nivel de rechazo estimado al 90%
  bbRT <- resRT[6]  # ... al 5%
  ccRT <- resRT[2]  # ... al 1%

     aaST=0
     bbST=0
     ccST=0
     aaRT=0
     bbRT=0
     ccRT=0

#  aTC <- (100-(which(resTC>0.10)[1]-1))/100
#  bTC <- (100-(which(resTC>0.05)[1]-1))/100
#  cTC <- (100-(which(resTC>0.01)[1]-1))/100
     
#  aaTC <- resTC[11] # Nivel de rechazo estimado al 90%
#  bbTC <- resTC[6]  # ... al 5%
#  ccTC <- resTC[2]  # ... al 1%
 
       
  vec <- correccion_NA(c(aRT,bRT,cRT))
  aRT <- vec[1]
  bRT <- vec[2] 
  cRT <- vec[3]

  vec <- correccion_NA(c(aST,bST,cST))
  aST <- vec[1]
  bST <- vec[2] 
  cST <- vec[3]

     
#  vec <- correccion_NA(c(aTC,bTC,cTC))
#  aTC <- vec[1]
#  bTC <- vec[2] 
#  cTC <- vec[3]
     
 
 
  # Almacenamiento de resultados
#  RESULTADOS$TC.90a[i]=aTC
#  RESULTADOS$TC.95a[i]=bTC
#  RESULTADOS$TC.99a[i]=cTC
#  RESULTADOS$TC.90b[i]=aaTC
#  RESULTADOS$TC.95b[i]=bbTC
#  RESULTADOS$TC.99b[i]=ccTC
  
  RESULTADOS$RT.90a[i]=aRT
  RESULTADOS$RT.95a[i]=bRT
  RESULTADOS$RT.99a[i]=cRT
  RESULTADOS$RT.90b[i]=aaRT
  RESULTADOS$RT.95b[i]=bbRT
  RESULTADOS$RT.99b[i]=ccRT

  RESULTADOS$ST.90a[i]=aST
  RESULTADOS$ST.95a[i]=bST
  RESULTADOS$ST.99a[i]=cST
  RESULTADOS$ST.90b[i]=aaST
  RESULTADOS$ST.95b[i]=bbST
  RESULTADOS$ST.99b[i]=ccST 



     
  # Impresion de resultados
  if(impresion){
     x11()
     par(mfrow=c(2,2))
 
  plot(resAD,main=paste(generador,"- Anderson Darling n =",n),
       sub="100 muestras de n datos", ylab="A*A",xlab="")
    lines(c(60,100),c(1.933,1.933),col=3,lwd=3)
    lines(c(60,100),c(2.492,2.492),col=2,lwd=3)
    lines(c(60,100),c(3.857,3.857),col=4,lwd=3)  
    legend("topleft",lty=1,lwd=3,col=c(3,2,4),
           legend=c( paste( aAD , "vs 0.90"),
                     paste( bAD , "vs 0.95"),
                     paste( cAD , "vs 0.99")
           )
    )
 
  plot(resKS,main=paste(generador,"- Kolmogorov Smirnov n =",n),
       sub="100 muestras de n datos", ylab="",xlab="")
    lines(c(0,20),rep(0.10,2),col=3,lwd=3)
    lines(c(0,20),rep(0.05,2),col=2,lwd=3)
    lines(c(0,20),rep(0.01,2),col=4,lwd=3)
 
     legend("bottomright",
           legend=c( paste( round(aaKS,5),",", round(aKS,3), "vs 0.90"),
                     paste( round(bbKS,5),",", round(bKS,3), "vs 0.95"),
                     paste( round(ccKS,5),",", round(cKS,3), "vs 0.99")
                   )
           )
 
  plot(resX2,main=paste(generador,"- X2 de Pearson n =",n),
       sub="100 muestras de n datos", ylab="",xlab="")
      lines(c(0,20),rep(0.10,2),col=3,lwd=3)
      lines(c(0,20),rep(0.05,2),col=2,lwd=3)
      lines(c(0,20),rep(0.01,2),col=4,lwd=3)
      legend("bottomright",
           legend=c( paste( round(aaX2,5),",", round(aX2,3), "vs 0.90"),
                     paste( round(bbX2,5),",", round(bX2,3), "vs 0.95"),
                     paste( round(ccX2,5),",", round(cX2,3), "vs 0.99")
           )
      )

  plot(sort(resMC1),main=paste(generador,"- Monte Carlo n =",n),
       sub="100 muestras de n datos", ylab="",xlab="",type="n",xlim=c(0,200),ylim=c(0,1.2))
      lines(c(0,20),rep(1/3,2),col=3,lwd=3)
      lines(c(0,20),rep(0.5,2),col=2,lwd=3)
      lines(c(0,20),rep(1,2),col=4,lwd=3)
      points(sort(resMC1),col=3,lwd=3)
      points(sort(resMC2),col=2,lwd=3)
      points(sort(resMC3),col=4,lwd=3)
      legend("right",
           legend=c( paste( "Funcion:",round(mean(resMC1),4),"vs 1/3"),
                     paste( "  Var=", round(var(resMC1),3)),
                     paste( "Tablero:",round(mean(resMC2),5),"vs 0.5"),
                     paste( "  Var=", round(var(resMC2),3)),
                     "Conteos Poisson"," por cuadro negro:",
                     paste( "Med=", round(mean(resMC3),3)*100,"%"),
                     paste( "  Var=", round(var(resMC3),3)*100, "%" )
           )
      )

 
        x11()
     par(mfrow=c(2,1))
 
 
 
  plot(resRT,main=paste(generador,"- Test de Rachas n =",n),
       sub="100 muestras de n datos", ylab="",xlab="")
      lines(c(0,20),rep(0.10,2),col=3,lwd=3)
      lines(c(0,20),rep(0.05,2),col=2,lwd=3)
      lines(c(0,20),rep(0.01,2),col=4,lwd=3)
      legend("bottomright",
           legend=c( paste( round(aaRT,5),",", round(aRT,3), "vs 0.90"),
                     paste( round(bbRT,5),",", round(bRT,3), "vs 0.95"),
                     paste( round(ccRT,5),",", round(cRT,3), "vs 0.99")
           )
      )


  plot(resST,main=paste(generador,"- Test Serial n =",n),
       sub="100 muestras de n datos", ylab="",xlab="")
      lines(c(0,20),rep(0.10,2),col=3,lwd=3)
      lines(c(0,20),rep(0.05,2),col=2,lwd=3)
      lines(c(0,20),rep(0.01,2),col=4,lwd=3)
      legend("bottomright",
          legend=c( paste( round(aaST,5),",", round(aST,3), "vs 0.90"),
                     paste( round(bbST,5),",", round(bST,3), "vs 0.95"),
                     paste( round(ccST,5),",", round(cST,3), "vs 0.99")
           )
      )


 
     
   }
  
   }
    
    
     # Lo siguiente se prueba con 100 muestras de 1000,
     #   con 10 muestras de 10 000
     #   y 1 muestra de 100 000
    
    for(n in c(1000,2000,5000,10000) ) { 
      ## APLICACION DE LOS TEST DE INDEPENDENCIA
     
      if(n==1000)
         for(i in 1:100){
           mm <- muestra[((i-1)*n+1):(n*i)]
           resRT[i]=Prueba.Rachas(mm)
           resST[i]=Prueba.Serial(mm)
#           resTC[i]=Tablas.Contingencia(mm)
         }
 
      if(n==10000){
         for(i in 1:10){
           mm <- muestra[((i-1)*n+1):(n*i)]
           resRT[i]=Prueba.Rachas(mm)
           resST[i]=Prueba.Serial(mm)
#           resTC[i]=Tablas.Contingencia(mm)
         }
         resRT  <- resRT[1:10]
         resST  <- resST[1:10]

#         resTC  <- resTC[1]
      }
      if(n==20000){
         for(i in 1:10){
           mm <- muestra[((i-1)*n+1):(n*i)]
           resRT[i]=Prueba.Rachas(mm)
           resST[i]=Prueba.Serial(mm)
#           resTC[i]=Tablas.Contingencia(mm)
         }
         resRT  <- resRT[1:5]
         resST  <- resST[1:5]

#         resTC  <- resTC[1]
      }

      if(n==100000){
         mm <- muestra
         #resRT=Prueba.Rachas(mm) # Regresa NA
         resRT=0
         resST=Prueba.Serial(mm)
#        resTC=Tablas.Contingencia(mm)
      }
   
    
 
     resRT=sort(resRT)
     resST=sort(resST)
#     resTC=sort(resTC)
 
 
  aRT <- (100-(which(resRT>0.10)[1]-1))/100
  bRT <- (100-(which(resRT>0.05)[1]-1))/100
  cRT <- (100-(which(resRT>0.01)[1]-1))/100
 
  aST <- (100-(which(resST>0.10)[1]-1))/100
  bST <- (100-(which(resST>0.05)[1]-1))/100
  cST <- (100-(which(resST>0.01)[1]-1))/100
 
  if(n!=100000){
  aaST <- resST[11] # Nivel de rechazo estimado al 90%
  bbST <- resST[6]  # ... al 5%
  ccST <- resST[2]  # ... al 1%

  aaRT <- resRT[11] # Nivel de rechazo estimado al 90%
  bbRT <- resRT[6]  # ... al 5%
  ccRT <- resRT[2]  # ... al 1%
  }
  if(n==100000){
     aaST=0
     bbST=0
     ccST=0
     aaRT=0
     bbRT=0
     ccRT=0
  }

#  aTC <- (100-(which(resTC>0.10)[1]-1))/100
#  bTC <- (100-(which(resTC>0.05)[1]-1))/100
#  cTC <- (100-(which(resTC>0.01)[1]-1))/100
     
#  aaTC <- resTC[11] # Nivel de rechazo estimado al 90%
#  bbTC <- resTC[6]  # ... al 5%
#  ccTC <- resTC[2]  # ... al 1%
 
       
  vec <- correccion_NA(c(aRT,bRT,cRT))
  aRT <- vec[1]
  bRT <- vec[2] 
  cRT <- vec[3]

  vec <- correccion_NA(c(aST,bST,cST))
  aST <- vec[1]
  bST <- vec[2] 
  cST <- vec[3]

  vec <- correccion_NA(c(aRT,bRT,cRT))
  aaRT <- vec[1]
  bbRT <- vec[2] 
  ccRT <- vec[3]

     
#  vec <- correccion_NA(c(aTC,bTC,cTC))
#  aTC <- vec[1]
#  bTC <- vec[2] 
#  cTC <- vec[3]
     
 
 
  # Almacenamiento de resultados
  i <- which(c(1000,5000)==n)[1]
#  RESULTADOS_INDEP$TC.90a[i]=aTC
#  RESULTADOS_INDEP$TC.95a[i]=bTC
#  RESULTADOS_INDEP$TC.99a[i]=cTC
#  RESULTADOS_INDEP$TC.90b[i]=aaTC
#  RESULTADOS_INDEP$TC.95b[i]=bbTC
#  RESULTADOS_INDEP$TC.99b[i]=ccTC
  
  RESULTADOS_INDEP$RT.90a[i]=aRT
  RESULTADOS_INDEP$RT.95a[i]=bRT
  RESULTADOS_INDEP$RT.99a[i]=cRT
  RESULTADOS_INDEP$RT.90b[i]=aaRT
  RESULTADOS_INDEP$RT.95b[i]=bbRT
  RESULTADOS_INDEP$RT.99b[i]=ccRT

  RESULTADOS_INDEP$ST.90a[i]=aST
  RESULTADOS_INDEP$ST.95a[i]=bST
  RESULTADOS_INDEP$ST.99a[i]=cST
  RESULTADOS_INDEP$ST.90b[i]=aaST
  RESULTADOS_INDEP$ST.95b[i]=bbST
  RESULTADOS_INDEP$ST.99b[i]=ccST 

  write.csv(RESULTADOS_INDEP,"Res2temp.csv")
 
 
  # Impresion de resultados
  if(impresion){
     x11()
     par(mfrow=c(2,1))
 


  plot(resST,main=paste(generador,"- Test Serial n =",n),
       sub="100 muestras de n datos", ylab="",xlab="")
      lines(c(0,20),rep(0.10,2),col=3,lwd=3)
      lines(c(0,20),rep(0.05,2),col=2,lwd=3)
      lines(c(0,20),rep(0.01,2),col=4,lwd=3)
      legend("bottomright",
          legend=c( paste( round(aaST,5),",", round(aST,3), "vs 0.90"),
                     paste( round(bbST,5),",", round(bST,3), "vs 0.95"),
                     paste( round(ccST,5),",", round(cST,3), "vs 0.99")
           )
      )
 
 
  plot(resRT,main=paste(generador,"- Test de Rachas n =",n),
       sub="100 muestras de n datos", ylab="",xlab="")
      lines(c(0,20),rep(0.10,2),col=3,lwd=3)
      lines(c(0,20),rep(0.05,2),col=2,lwd=3)
      lines(c(0,20),rep(0.01,2),col=4,lwd=3)
      legend("bottomright",
           legend=c( paste( round(aaRT,5),",", round(aRT,3), "vs 0.90"),
                     paste( round(bbRT,5),",", round(bRT,3), "vs 0.95"),
                     paste( round(ccRT,5),",", round(cRT,3), "vs 0.99")
           )
      )




#  plot(resTC,main=paste(generador,"- T. Contingencia dim=2, n =",n),
#       sub="100 muestras de n datos", ylab="",xlab="")
#    lines(c(0,20),rep(0.10,2),col=3,lwd=3)
#    lines(c(0,20),rep(0.05,2),col=2,lwd=3)
#    lines(c(0,20),rep(0.01,2),col=4,lwd=3)
 
#     legend("bottomright",
#           legend=c( paste( round(aaTC,5),",", round(aTC,3), "vs 0.90"),
#                     paste( round(bbTC,5),",", round(bTC,3), "vs 0.95"),
#                     paste( round(ccTC,5),",", round(cTC,3), "vs 0.99")
#                   )
#           )
 
     
   }
   }
   return(RESULTADOS)
} 
 
  correccion_NA <- function(vec){
     if(is.na(vec[3])){
       vec[3]=1
     }
     if(is.na(vec[2])){
        vec[3]=1
        vec[1]=1
     }
     if(is.na(vec[1])){
        vec[3]=1
        vec[2]=1
        vec[1]=1
     }
 
    
     return(vec)
  }   
 
 
 
#############################################  Tablas de contingencia
 
Tabla_Cont<-function(muestra, particion=2*floor((sqrt(length(muestra))/2))) {
   k<-floor(length(muestra)/2)
   A <- matrix(0,particion,particion)
   R<-particion*muestra
   for(i in 1:k){
      a<-floor(R[i,1])
      b<-floor(R[i,2])
      A[a+1,b+1]<-A[a+1,b+1]+1
   }  
   B<-((A-((k)/(particion^2)))^2)/((k)/(particion^2))
   CHI<-sum(B)
   return(1-pchisq(CHI,df=(particion-1)*(particion-1)))
}
  
Tablas.Contingencia<-function(muestra,particion=5,impresion=F) {
 
   k<-floor(length(muestra)/2)
   Muestra<-cbind(muestra[1:k],muestra[(k+1):(2*k)])
 
 
   X<-c(1:particion)
   k<-floor(length(Muestra)/(2*particion))
   for(i in 1:particion)
   {
      X[i]<-Tabla_Cont(Muestra[((i-1)*k):(i*k),])
   }
   if(impresion){
      plot(1:particion,X)}
   return(sum(X)/particion)
}




   ##################################################
   # Mediciones con los resultados
 
   # Se toma el promedio de los resultados
   #    para los distintos tamaños de muestra
   #    y se calcula un intervalo con +/- una vez
   #    la siguiente medida de variación:
   # Si se rechazan más pruebas que las
   #    esperadas, se le da el doble del peso
   #    a que si no se rechazan tantas pruebas
   #    como las esperadas a rechazar.
 
  
   Medida <- function(vector,valor_esperado,arriba=T){
     # Arriba=T es que no se rechaza para arriba
       media <- mean(vector)
      
       v <- vector-valor_esperado
 
       if(arriba){
         positivos <- as.numeric(v<0)+1
       }
       if(!arriba){
         positivos <- as.numeric(v>0)+1
       }
      
       v  <- (v^2)*positivos/3
       v  <- sum(v)/length(v)
       variacion <- v[1]  # Variación
      
       return( c(media,variacion) )
   }
  
   Medicion <- function(R1){
  
      R2 <- R1[1,]
      R2[2,] <- R1[1,]
     
      # Valores de Monte Carlo
        # Para la variacion: variacion con lo esperado
        #   como media ponderada junto con las variaciones
        #   para n=1000, 10 000
     
        y <- na.omit(R1[1:2,21])
        R2[1,20]=mean(na.omit(R1[1:2,20]))
        R2[2,20]=(sum(y)/length(y)+2*sum((R1[1:3,20]-1/3)^2)/3 )/3
        R2[1,22]=mean(R1[1:2,22])
        y <- na.omit(R1[1:2,23])
        R2[2,22]=(sum(y)/length(y)+2*sum((R1[1:3,22]-1/3)^2)/3 )/3
        # La variacion está ponderada por 1/3(var estadistcos)
        #   y 2/3 (var resultados con respecto a lo esperado)
     
      R2[,1] <- c("Media","Variacion")
      names(R2)[1] <- ""
 
 
      R2[,5] <- Medida(R1[,5],1.933,F)
      R2[,6] <- Medida(R1[,6],2.492,F)
      R2[,7] <- Medida(R1[,7],3.857,F)
 
 
      for(i in c(2,8,14)){
         R2[,i]   <- Medida(R1[,i],  0.90,T)
         R2[,i+1] <- Medida(R1[,i+1],0.95,T)
         R2[,i+2] <- Medida(R1[,i+2],0.99,T)
      }
 
      for(i in c(11,17)){
        R2[,i]   <- Medida(R1[,i],  0.90,F)
        R2[,i+1] <- Medida(R1[,i+1],0.95,F)
        R2[,i+2] <- Medida(R1[,i+2],0.99,F)
      }

      i=20
        R2[,i]   <- Medida(R1[,i],  1/3,F)
        R2[,i+1] <- mean(R1[,i+1])
        R2[,i+2] <- Medida(R1[,i+2],  0.5,F)
        R2[,i+3] <- mean(R1[,i+3])
        R2[,i+4] <- Medida(R1[,i+4],  1,F)
        R2[,i+5] <- mean(R1[,i+5])

      for(i in 26:36){
        R2[1,i] <- mean(R1[,i])
        R2[2,i] <- var(R1[,i])
      }



       
      return(R2)
   }
 
 
   Valor_Generador_Unif <- function( R2 ){
 
   # Para dar un número que identifique a este generador
   #   se toman una media ponderada de las variaciones
   #   con peso 3/8 para AD, 2/8 KS, 1/8 X2 y 2/8 MC
      AA <- sum(R2[2,-1]*c( rep(3/8,6) , rep(2/8,6) , rep(1/8,6), rep(2/8,6) ))
      return (AA)
   }
 
   Valor_Generador_Indep <- function( R2 ){
 
   # Para dar un número que identifique a este generador
   #   se toman una media ponderada de las variaciones
   #   con peso 1/2 para RT, 1/2 ST    (y 0 TC, no mide indep)
 
   # Valor de Runs test con un tipo de varianza ponderada
      var_RT <- c((R2$RT.90a[1]-0.90), (R2$RT.95a[1]-0.95), (R2$RT.99a[1]-0.99))
      var_RT <- sum( (var_RT^2)/c(0.9,0.95,0.99) )*2/3
      var_RT <- var_RT+sum(R2[2,2:7]*1/3)
 
 
   # Valor de Serial test con un tipo de varianza ponderada
      var_ST <- c((R2$ST.90a[1]-0.90), (R2$ST.95a[1]-0.95), (R2$ST.99a[1]-0.99))
      var_ST <- sum( (var_ST^2)/c(0.9,0.95,0.99) )*2/3
      var_ST <- var_ST+sum(R2[2,8:13]*1/3)
 
 
   # Valor de las tablas de contingencia
   #   var_TC <- (((R2$TC.90a[1]-0.9)^2)/0.9+((R2$TC.95a[1]-0.95)^2)/0.95+((R2$TC.99a[1]-0.99)^2)/0.99 )/3
   #   var_TC <- (2*var_TC+sum(c(R2$TC.90a[2],R2$TC.95a[2],R2$TC.99a[2] ))/3)/3
 
      return ( var_RT*0.5+var_ST*0.5 )
   }
  
 
   filtro_na <- function( tabla ){
      nr <- dim(tabla)[1]
      nc <- dim(tabla)[2]
     
      for(i in 1:nc){
         tabla[which(is.na(tabla[,i])),i]=0
      }
      return(tabla)
   }



########################################################
##     Compilar primero de aqui hasta el inicio        #
########################################################
##       despues se puede compilar paso a paso para
#        explorar los valores de los objetos


#
##
########################################################
##                      EJECUCION                      #
########################################################
 
  
########################################################
########################################   Base de datos
 
 
# Tamaños de las muestras
   tam <- c(20,40,100,200,500,1000)
  
# 100,000 aleatorios con misma semilla; de R, C y H.
# 100,000 aleatorios con distinta semilla cada mil (R y C)
  
   set.seed(33)
   Muestra_R  <- runif(100000,0,1)
   Muestra_R2 <- c()
  
   for(i in 1:100){
     set.seed(i*i)
     Muestra_R2 <- c(Muestra_R2,runif(1000,0,1))
   }
  
   write.table( Muestra_R, "Aleat_R.txt" )
   write.table( Muestra_R2,"Aleat_R2.txt" )
 
  
   Muestra_C  <- scan("Aleat_C.txt")  # Semilla=33
   Muestra_C2 <- scan("Aleat_C2.txt") # Varia la semilla cada mil y es 0,1000,2000,...
   Muestra_H  <- read.table("Aleat_H.txt")[,1]
  
      if(!is.numeric(Muestra_C)){
         Muestra_C <- as.numeric(Muestra_C)
      }
      if(!is.numeric(Muestra_C2)){
         Muestra_C2 <- as.numeric(Muestra_C2)
      }
      if(!is.numeric(Muestra_H)){
         Muestra_H <- as.numeric(Muestra_H)
      }
 
   Muestra_C  <- ((Muestra_C)%%99+1)/100 #Forma tradicional
   Muestra_C2 <- ((Muestra_C2)%%99+1)/100

     # Dividiendo así a C, se preservan más dígitos de los números

   Muestra_H  <- as.vector( (Muestra_H+1)/257 )
  
   # Eran enteros del 0 al 255
  
   # Ya se tienen las muestras con valores en (0,1)
 
  
  
########################################################
#################################   Pruebas y Mediciones
     
  
   R1_res  <- Prueba(Muestra_R, tam,"R",T)
   R1_res2 <- read.csv("Res2temp.csv")
   R2_res  <- Prueba(Muestra_R2,tam,"R2",T)
   R2_res2 <- read.csv("Res2temp.csv")
   C1_res  <- Prueba(Muestra_C, tam,"C",T)
   C1_res2 <- read.csv("Res2temp.csv")
   C2_res  <- Prueba(Muestra_C2,tam,"C2",T)
   C2_res2 <- read.csv("Res2temp.csv")
   H1_res  <- Prueba(Muestra_H, tam,"H",T)
   H1_res2 <- read.csv("Res2temp.csv")
 
   par(mfrow=c(1,1))
  
   # Guardado en archivo
   write.csv(R1_res,"Resultados_R1.csv")
   write.csv(R2_res,"Resultados_R2.csv")
   write.csv(C1_res,"Resultados_C1.csv")
   write.csv(C2_res,"Resultados_C2.csv")
   write.csv(H1_res,"Resultados_H1.csv")

   write.csv(R1_res2,"Resultados2_R1.csv")
   write.csv(R2_res2,"Resultados2_R2.csv")
   write.csv(C1_res2,"Resultados2_C1.csv")
   write.csv(C2_res2,"Resultados2_C2.csv")
   write.csv(H1_res2,"Resultados2_H1.csv")

   R1_res <- filtro_na(R1_res)
   R2_res <- filtro_na(R2_res)
   C1_res <- filtro_na(C1_res)
   C2_res <- filtro_na(C2_res)
   H1_res <- filtro_na(H1_res)

   R1_res2 <- filtro_na(R1_res2)
   R2_res2 <- filtro_na(R2_res2)
   C1_res2 <- filtro_na(C1_res2)
   C2_res2 <- filtro_na(C2_res2)
   H1_res2 <- filtro_na(H1_res2)


   R1_med <- Medicion(R1_res)
   R2_med <- Medicion(R2_res)
   C1_med <- Medicion(C1_res)
   C2_med <- Medicion(C2_res)
   H1_med <- Medicion(H1_res)
 
########################################################
################################   Resultados Prefinales
prefinales <- function(R1,tit=""){
   x11()
   par(mfrow=c(2,2))
   plot( c(0,1),c(0.9,0.9), type="l",lwd=2,main=paste(tit,"Pruebas al 0.90" ))
      points(seq(0,1,length=5),as.vector(R1[1,c(3,9,15,25,31)]), col=c(2,3,4,7), lwd=5)
      legend("topright", col=c(2,3,4,7), lwd=3,legend=c("AD","KS","X2","RT","TC") )
   plot( c(0,1),c(0.95,0.95), type="l",lwd=2,main=paste(tit,"Pruebas al 0.95") )
      points(seq(0,1,length=5),as.vector(R1[1,c(4,10,16,26,32)]), col=c(2,3,4,7), lwd=5)
      legend("topright", col=c(2,3,4,7), lwd=3,legend=c("AD","KS","X2","RT","TC") )
   plot( c(0,1),c(0.99,0.99), type="l",lwd=2,main=paste(tit,"Pruebas al 0.99" ))
      points(seq(0,1,length=5),as.vector(R1[1,c(5,11,17,27,33)]), col=c(2,3,4,7), lwd=5)
      legend("topright", col=c(2,3,4,7), lwd=3,legend=c("AD","KS","X2","RT","TC") )
   plot( c(0,1),c(0.99,0.99), type="n",lwd=2,main=paste(tit,"Monte Carlo" ))
      legend("right", legend=c(paste(round(R1[1,21],4),"vs",0.3333),
                               paste(" variacion: ",round(R1[1,22],6)),
                               paste(round(R1[1,23],4),"vs",0.5000),
                               paste(" variacion: ",round(R1[1,24],6))
             ))
}
 
 
   prefinales(R1_med,"R1")
   prefinales(R2_med,"R2")
   prefinales(C1_med,"C1")
   prefinales(C2_med,"C2")
   prefinales(H1_med,"H1")


 
########################################################
######################################   Resultado Final
 
 
   R1_val <- Valor_Generador_Unif(R1_med)
   R2_val <- Valor_Generador_Unif(R2_med)
   C1_val <- Valor_Generador_Unif(C1_med)
   C2_val <- Valor_Generador_Unif(C2_med)
   H1_val <- Valor_Generador_Unif(H1_med)
 
   R1_val2 <- Valor_Generador_Indep(R1_med)
   R2_val2 <- Valor_Generador_Indep(R2_med)
   C1_val2 <- Valor_Generador_Indep(C1_med)
   C2_val2 <- Valor_Generador_Indep(C2_med)
   H1_val2 <- Valor_Generador_Indep(H1_med)
 
 
   ###########################
   # Resultados de uniformidad
  
   x11()
   par(mfrow=c(1,2))
   Res <- c(R1_val,R2_val,C1_val,C2_val,H1_val)
   plot(as.factor(c("R1","R2","C1","C2","H1")),Res,
        main="Valor de uniformidad")
   legend("right",legend=c( paste("R1: ",round(Res[1],3)),paste("R2: ",round(Res[2],3)),paste("H1: ",round(Res[5],3)),paste("C1: ",round(Res[3],3)),paste("C2: ",round(Res[4],3)) ))
 
   # Para dar una medicion de la confianza en
   #   cada generador
   Res <- 1-Res/sum(Res)
   plot(as.factor(c("R1","R2","C1","C2","H1")),Res,
        main="Confianza de los generadores")
   legend("right",legend=c( paste("R1: ",round(Res[1],3)),paste("R2: ",round(Res[2],3)),paste("H1: ",round(Res[5],3)),paste("C1: ",round(Res[3],3)),paste("C2: ",round(Res[4],3)) ))
 
 
   ############################
   # Resultado de independencia
 
   x11()
   par(mfrow=c(1,2))
   Res <- c(R1_val2,R2_val2,C1_val2,C2_val2,H1_val2)
   plot(as.factor(c("R1","R2","C1","C2","H1")),Res,
        main="Valor de independencia")
   legend("right",legend=c( paste("R1: ",round(Res[1],3)),paste("R2: ",round(Res[2],3)),paste("H1: ",round(Res[5],3)),paste("C1: ",round(Res[3],3)),paste("C2: ",round(Res[4],3)) ))


########################################################
##########################   Confianza en cada generador
 
   # Para dar una medicion de la confianza en cada generador
   Res <- 1-Res/sum(Res)
   plot(as.factor(c("R1","R2","C1","C2","H1")),Res,
        main="Confianza de los generadores")
   legend("right",legend=c( paste("R1: ",round(Res[1],3)),paste("R2: ",round(Res[2],3)),paste("H1: ",round(Res[5],3)),paste("C1: ",round(Res[3],3)),paste("C2: ",round(Res[4],3)) ))
 
 
  
 
# Fin
