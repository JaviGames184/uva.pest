#### funciones auxiliares

clave<-c(toupper(letters)[1:14],"Ñ",toupper(letters)[15:26]," ")

letra_a_numero<-function(a) match(a,clave) 
numero_a_letra<-function(n) clave[n]

## funcion que transforma mensaje de texto 
## a secuencia de numeros

mensaje_a_numero<-function(mensaje){
  mensajenumerado<-NULL
  for(i in 1:nchar(mensaje)){
    mensajenumerado<-c(mensajenumerado,letra_a_numero(substring(mensaje,i,i)))  
  }
  return(mensajenumerado)
}

## funcion que transforma secuencia de numeros
## a mensaje de texto 

mensaje_a_letra<-function(mensaje){
  mensajetexto<-NULL
  for(i in 1:length(mensaje)){
    mensajetexto<-paste(mensajetexto,numero_a_letra(mensaje[i]),sep='')  
  }
  return(mensajetexto)
}

bigramas_palabra <- function(palabra) {
  aaa<-unlist(strsplit(palabra,''))
  nnn<-length(aaa)
  bigramas<-paste0(aaa[1:(nnn-1)],aaa[2:nnn])
  return(bigramas)
}

trigramas_palabra <- function(palabra) {
  aaa<-unlist(strsplit(palabra,''))
  nnn<-length(aaa)
  trigramas<-paste0(aaa[1:(nnn-2)],aaa[2:(nnn-1)],aaa[3:nnn])
  return(trigramas)
}



## funcion que aplica el descifrado de sustitucion 
## dado por ordenacion

decodificador<-function(mensaje,ordenacion){
  fun.auxiliar<-function(i) ordenacion[i]
  mensaje_decodificado<-NULL
  for (j in 1:length(mensaje)) mensaje_decodificado<-c(mensaje_decodificado,fun.auxiliar(mensaje[j]))
  return(mensaje_decodificado)
}

