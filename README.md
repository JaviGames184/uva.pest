# Descifrar un mensaje
Esta práctica tiene como objetivo utilizar el algoritmo Metropolis-Hastings para descifrar el mensaje:

HTSVBPX QTIVBHRXIXBSX XLNTBVBGVBGREX HVSBSXBPXIQVYRXIHTBSXBLTILRXILRVBCBSXB XGRJRTIBXQHXBSX XLNTBRILGACXBGVBGREX HVSBSXBLVYERV BSXB XGRJRTIBTBSXBL XXILRVBVQRBLTYTBGVBGREX HVSBSXBYVIRWXQHV BQAB XGRJRTIBTBQABL XXILRVBRISRZRSAVGBCBLTGXLHRZVYXIHXBHVIHTBXIBPAEGRLTBLTYTBXIBP RZVSTBPT BGVBXIQXMVIDVBGVBP VLHRLVBXGBLAGHTBCBGVBTEQX ZVILRVBHTSTBRISRZRSATBHRXIXBSX XLNTBVBGVBGREX HVSBSXBTPRIRTIBCBSXBXKP XQRTIBXQHXBSX XLNTBRILGACXBXGBSXBITBQX BYTGXQHVSTBVBLVAQVBSXBQAQBTPRIRTIXQBXGBSXBRIZXQHRJV BCB XLRER BRIWT YVLRTIXQBCBTPRIRTIXQBCBXGBSXBSRWAISR GVQBQRIBGRYRHVLRTIBSXBW TIHX VQBPT BLAVGUARX BYXSRTBSXBXKP XQRTI

# Solución propuesta
## Ficheros de datos

En primer lugar importamos los diferentes archivos con la información que vamos a necesitar:

-   bigramas_espanol: Archivo con las frecuencias correspondientes a los bigramas

-   trigramas_espanol: Archivo con las frecuencias correspondientes a los trigramas

-   long_palabras_espanol: Archivo con las frecuencias correspondientes a las longitudes de las palabras en español

```{r}
bigramas_espanol <- read.csv("D:/Archivos/4º/1Q/PEST/Entregable1/DecodificacionMCMC/bigramas_espanol.txt", sep="", encoding="UTF-8")
# Para cambiar el bigrama "NA":
bigramas_espanol[28,1] <- "NA"

trigramas_espanol <- read.csv("D:/Archivos/4º/1Q/PEST/Entregable1/DecodificacionMCMC/trigramas_espanol.txt", sep="", encoding="UTF-8")

long_palabras_espanol <- read.csv("D:/Archivos/4º/1Q/PEST/Entregable1/DecodificacionMCMC/long_palabras_espanol.txt", sep="", encoding="UTF-8")
```

## Funciones auxiliares

Introducimos también las funciones que se nos proporcionan para así poder utilizarlas más adelante

```{r}
clave<-c(toupper(letters)[1:14],"Ñ",toupper(letters)[15:26]," ")
letra_a_numero<-function(a) match(a,clave) 
numero_a_letra<-function(n) clave[n]
```

-   mensaje_a_numero: Función que transforma un mensaje de texto a una secuencia de números

    ```{r}
    mensaje_a_numero<-function(mensaje){
      mensajenumerado<-NULL
      for(i in 1:nchar(mensaje)){
        mensajenumerado<-c(mensajenumerado,letra_a_numero(substring(mensaje,i,i)))  
      }
      return(mensajenumerado)
    }
    ```

-   mensaje_a_letra: Función que transforma una secuencia de números a un mensaje de texto

    ```{r}
    mensaje_a_letra<-function(mensaje){
      mensajetexto<-NULL
      for(i in 1:length(mensaje)){
        mensajetexto<-paste(mensajetexto,numero_a_letra(mensaje[i]),sep='')  
      }
      return(mensajetexto)
    }
    ```

-   bigramas_palabra: Función que extrae los bigramas de una palabra

    ```{r}
    bigramas_palabra <- function(palabra) {
      aaa<-unlist(strsplit(palabra,''))
      nnn<-length(aaa)
      bigramas<-paste0(aaa[1:(nnn-1)],aaa[2:nnn])
      return(bigramas)
    }
    ```

-   trigramas_palabra: Función que extrae los trigramas de una palabra

    ```{r}
    trigramas_palabra <- function(palabra) {
      aaa<-unlist(strsplit(palabra,''))
      nnn<-length(aaa)
      trigramas<-paste0(aaa[1:(nnn-2)],aaa[2:(nnn-1)],aaa[3:nnn])
      return(trigramas)
    }
    ```

-   decodificador: Función que aplica el descifrado de sustitucion dado por ordenacion

    ```{r}
    decodificador<-function(mensaje,ordenacion){
      fun.auxiliar<-function(i) ordenacion[i]
      mensaje_decodificado<-NULL
      for (j in 1:length(mensaje)) mensaje_decodificado<-c(mensaje_decodificado,fun.auxiliar(mensaje[j]))
      return(mensaje_decodificado)
    }
    ```


## Descifrado mediante el algoritmo Metropolis

### Creamos una función para calcular los Scores

Función para calcular el Score:

$$score_{j}(f) = \prod_{palabra \in f(mensaje)}s_{j}(palabra)$$

Siendo:

$score_{1}(palabra)=$ la frecuencia observada de la longitud de palabra

$score_{2}(palabra)=$ el producto de las frecuencias observadas de los bigramas en la palabra

$score_{3}(palabra)=$ el producto de las frecuencias observadas de los trigramas en la palabra

(En el caso de que uno de los valores no se encuentre en el archivo de bigramas o de trigramas se le asignará un valor arbitrario de $10^{-8}$)

A partir de los scores anteriores utilizamos $score(f) = (score_{1}(f))^{3}score_{2}(f)(score_{3}(f))^{0.5}$

<p> </p>

Como estos valores son muy grandes decidimos tomar logaritmos de forma que la fórmula que utilizamos es:

$$\ln(score_{j}(f)) = \sum_{palabra \in f(mensaje)} \ln(s_{j}(palabra))$$

Siendo:

$\ln(score_{1}(palabra))=$ el logaritmo de la frecuencia observada de la longitud de palabra

$\ln(score_{2}(palabra))=$ la suma del logaritmo de las frecuencias observadas de los bigramas en la palabra

$\ln(score_{3}(palabra))=$ la suma del logaritmo de las frecuencias observadas de los trigramas en la palabra

Al haber aplicado los logaritmos, el logaritmo del score total será: $$\ln(score(f)) = 3\times \ln(score_{1}(f))+\ln(score_{2}(f))+0.5\times \ln(score_{3}(f))$$

Creamos entonces la función que nos permite calcular el score anterior:

```{r Scores}
CalcularScore <- function(palabras){
palabras <- palabras[[1]]

# Score 1
longitud <- str_length(palabras)
frecuencias_encontradas <- long_palabras_espanol[long_palabras_espanol[,1] %in% longitud, ]
longitud_encontrada <- table(longitud)
longitud_encontrada <- longitud_encontrada[names(longitud_encontrada) %in% frecuencias_encontradas[,1]]
frecuencias <- rep(frecuencias_encontradas[,2], longitud_encontrada)
noencontrados <- length(longitud) - sum(longitud_encontrada)

score1 <- sum(log(frecuencias)) + (log(10^-8) * noencontrados)


# Score 2 
# Si la palabra es de una letra no tenemos bigramas
palabras_bigramas <- palabras[str_length(palabras) > 1]
bigramas <- bigramas_palabra(palabras_bigramas)
frecuencias_encontradas <- bigramas_espanol[bigramas_espanol[, 1] %in% bigramas, ]
# Ordenamos las frecuencias alfabéticamente
frecuencias_encontradas <- frecuencias_encontradas[order(frecuencias_encontradas[,1]), ]
bigrama_encontrado <- table(bigramas)
bigrama_encontrado <- bigrama_encontrado[names(bigrama_encontrado) %in% frecuencias_encontradas[,1]]
frecuencias <- rep(frecuencias_encontradas[,2], bigrama_encontrado)
noencontrados <- length(bigramas) - sum(bigrama_encontrado)
    
score2 <- sum(log(frecuencias)) + (log(10^-8) * noencontrados)


# Score 3
# Si la palabra es de una o dos letras no tenemos trigramas
palabras_trigramas <- palabras[str_length(palabras) > 2]
trigramas <- trigramas_palabra(palabras_trigramas)
frecuencias_encontradas <- trigramas_espanol[trigramas_espanol[, 1] %in% trigramas, ]
# Ordenamos las frecuencias alfabéticamente
frecuencias_encontradas <- frecuencias_encontradas[order(frecuencias_encontradas[,1]), ]
trigrama_encontrado <- table(trigramas)
trigrama_encontrado <- trigrama_encontrado[names(trigrama_encontrado) %in% frecuencias_encontradas[,1]]
frecuencias <- rep(frecuencias_encontradas[,2], trigrama_encontrado)
noencontrados <- length(trigramas) - sum(trigrama_encontrado)

score3 <- sum(log(frecuencias)) + (log(10^-8) * noencontrados)

# Calculamos el score total y lo devolvemos
score <- (3*score1) + score2 + (0.5*score3) 
return(score)
}
```

### El algoritmo Metropolis-Hastings

Introducimos el mensaje que queremos descifrar:

```{r mensaje_cifrado}
mensaje_cifrado <- "HTSVBPX QTIVBHRXIXBSX XLNTBVBGVBGREX HVSBSXBPXIQVYRXIHTBSXBLTILRXILRVBCBSXB XGRJRTIBXQHXBSX XLNTBRILGACXBGVBGREX HVSBSXBLVYERV BSXB XGRJRTIBTBSXBL XXILRVBVQRBLTYTBGVBGREX HVSBSXBYVIRWXQHV BQAB XGRJRTIBTBQABL XXILRVBRISRZRSAVGBCBLTGXLHRZVYXIHXBHVIHTBXIBPAEGRLTBLTYTBXIBP RZVSTBPT BGVBXIQXMVIDVBGVBP VLHRLVBXGBLAGHTBCBGVBTEQX ZVILRVBHTSTBRISRZRSATBHRXIXBSX XLNTBVBGVBGREX HVSBSXBTPRIRTIBCBSXBXKP XQRTIBXQHXBSX XLNTBRILGACXBXGBSXBITBQX BYTGXQHVSTBVBLVAQVBSXBQAQBTPRIRTIXQBXGBSXBRIZXQHRJV BCB XLRER BRIWT YVLRTIXQBCBTPRIRTIXQBCBXGBSXBSRWAISR GVQBQRIBGRYRHVLRTIBSXBW TIHX VQBPT BLAVGUARX BYXSRTBSXBXKP XQRTI"
```

Aplicamos el algoritmo Metropolis-Hastings de la siguiente forma:

Para comenzar el algoritmo necesitamos utilizar una configuración inicial que se escogerá de forma aleatoria.

A continuación, necesitamos generar un nuevo candidato j con probabilidad $q_{i,j}$. Para hacer esto vamos a elegir dos letras del mensaje anterior para intercambiarlas y de esta forma conseguir al nuevo candidato.

Por otra parte queremos que se acepte al nuevo candidato (nuevo mensaje donde hemos cambiado dos letras) con probabilidad $\alpha_{i,j}=min(1, \pi_jq_{j,i}/\pi_iq_{i,j})$. Nos podemos dar cuenta que esto en realidad sería la división del score calculado para el nuevo mensaje entre el score para el mensaje anterior dado que la matriz Q sería una matriz donde todos los valores serían iguales (salvo la donde $q_{i,i}=0$) pero como hemos tomado logaritmos esto sería: $$\alpha_{i,j}=min(1, e^{ln(scorenuevo)}/e^{ln(scoreantiguo)}) = min(1, e^{ln(scorenuevo) - ln(scoreantiguo)})$$.

La comprobación que se utilizará en el código será que si:

$min(1,e^{ln(scorenuevo) - ln(scoreantiguo)}) \le U(0,1) \Longrightarrow$ Nos quedamos con el mensaje antiguo

Creamos ahora un código que nos permita utilizar este algoritmo (siguiendo esos pasos):

```{r Metropolis-Hastings}
# Hacemos una reordenación aleatoria para empezar
nueva_ordenacion<-sample(1:28)
nuevo_mensaje <- mensaje_a_letra(decodificador(mensaje_a_numero(mensaje_cifrado), nueva_ordenacion))
# Calculamos el score para la configuración inicial
palabras <- strsplit(nuevo_mensaje, split=" ")
nuevo_score <- CalcularScore(palabras)

# Aplicamos el algoritmo de Metropolis-Hastings con 30000 iteraciones
start_time <- Sys.time()
n=30000
for (i in 1:n) {
	# Intercambiamos dos letras de forma aleatoria
    letras_cambiar <- sample(1:28,2)
    antigua_ordenacion <- nueva_ordenacion
    nueva_ordenacion[letras_cambiar[1]] <- antigua_ordenacion[letras_cambiar[2]]
    nueva_ordenacion[letras_cambiar[2]] <- antigua_ordenacion[letras_cambiar[1]]
    
    # Guardamos el nuevo mensaje y el anterior
    anterior_mensaje <- nuevo_mensaje
    nuevo_mensaje <- mensaje_a_letra(decodificador(mensaje_a_numero(mensaje_cifrado), nueva_ordenacion))
    # Calculamos el Score para el nuevo mensaje
    anterior_score <- nuevo_score
	palabras <- strsplit(nuevo_mensaje, split=" ")
	nuevo_score <- CalcularScore(palabras) 
	
    # Decidimos si aceptar al candidato generado
    if (min(exp(nuevo_score - anterior_score),1) <= runif(1)) {
        # Rechazamos el cambio de letras realizado y por tanto, volvemos al estado anterior
        nueva_ordenacion <- antigua_ordenacion
        nuevo_mensaje <- anterior_mensaje
        nuevo_score <- anterior_score
    }
}
end_time <- Sys.time()
```

```{r include=FALSE}
ordenacion_alt <- c(22,28,26,27,2,10,12,21,14,7,25,3,15,8,11,24,17,20,9,4,16,18,1,6,5,13,23,19)
nuevo_mensaje <- mensaje_a_letra(decodificador(mensaje_a_numero(mensaje_cifrado), ordenacion_alt))
palabras <- strsplit(nuevo_mensaje, split=" ")
nuevo_score <- CalcularScore(palabras)
```

Imprimimos los resultados que nos proporciona el algoritmo tras hacer todas las iteraciones y el tiempo que ha necesitado para hacerlo:

```{r Resultados}
# Duración de la ejecución
end_time-start_time

# Resultado final
mensaje_final <- nuevo_mensaje; mensaje_final
nuevo_score
```

Se ha necesitado hacer varias ejecuciones para poder llegar a este resultado ya que hemos visto que existen otras configuraciones que se pueden obtener como resultado además de la correcta.

Entre estas cabe destacar que una de ellas tiene un score aún mejor que el del mensaje correcto que consiste en hacer el intercambio de las letras "ñ" y "j". Tomando la ordenación que se tiene como resultado del algoritmo cuando este mensaje es el que se proporciona como respuesta vemos que este mensaje sería el siguiente:

```{r Resultado_Alternativo}
ordenacion_alt <- c(22,28,26,27,2,15,12,21,14,7,25,3,10,8,11,24,17,20,9,4,16,18,1,6,5,13,23,19)

# Resultado alternativo
mensaje_alt <- mensaje_a_letra(decodificador(mensaje_a_numero(mensaje_cifrado), ordenacion_alt))
mensaje_alt
palabras <- strsplit(mensaje_alt, split=" ")
score_alt <- CalcularScore(palabras)
score_alt
```

Vemos que el score de este mensaje alternativo es ligeramente superior al anterior pese a no ser el mensaje correcto (13737.83 - 13736.6 = 1.23) y por ese motivo en algunas de las ejecuciones del algoritmo es el que aparece como solución final.

