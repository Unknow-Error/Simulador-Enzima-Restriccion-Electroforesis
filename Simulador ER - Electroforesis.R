library(stringr)
library(glue)
library(Biostrings)

# Script realizado por Nicolás Salvatore - UNQ.
# Script que simula corte de enzima de restricción, en particular, la enzima hinfi y gráfica su electroforesis en gel de agarosa.
# Librerías requeridas: stringr, glue y Biostrings
# secuencia_adn: una cadena de texto que representa la secuencia de ADN
# sitio_restriccion: una cadena de texto que representa el sitio de corte de la enzima de restricción
# corte_pos: la posición dentro del sitio de restricción donde la enzima corta (1-indexed)

simular_enzima_restriccion <- function(secuencia_adn, sitio_restriccion, corte_pos) {
  # Inicializar una lista para almacenar los fragmentos de ADN
  fragmentos <- list()
  
  # Verificar que el sitio de restricción está contenido en la secuencia de ADN, si no lo está, retorna su secuencia sin cortar
  if (grepl(sitio_restriccion, secuencia_adn)) {
    # Encontrar todas las posiciones del sitio de restricción en la secuencia de ADN
    posiciones <- gregexpr(sitio_restriccion, secuencia_adn)[[1]]
    
    # Inicializar una lista para almacenar los fragmentos de ADN
    fragmentos <- list()
    
    # Variable para rastrear la posición inicial del siguiente fragmento
    inicio_fragmento <- 1
    
    # Recorrer cada posición donde se encuentra el sitio de restricción
    for (pos in posiciones) {
      # La posición de corte en la secuencia de ADN
      corte_real <- pos + corte_pos - 1
      
      # Extraer el fragmento de ADN hasta el punto de corte
      fragmento <- substr(secuencia_adn, inicio_fragmento, corte_real)
      fragmentos <- c(fragmentos, fragmento)
      
      # Actualizar la posición inicial del siguiente fragmento
      inicio_fragmento <- corte_real + 1
    }
    
    # Agregar el fragmento final desde el último punto de corte hasta el final de la secuencia
    fragmentos <- c(fragmentos, substr(secuencia_adn, inicio_fragmento, nchar(secuencia_adn)))  
  }
  else{
     fragmentos <- list(secuencia_adn)
  }
  return(fragmentos)
}

# Función para una enzima restricción con "N" en su secuencia con hinfi:
simular_corte_hinfi <- function(secuencia_adn, corte_pos) {
  nucleotidos <- c("A", "T", "C", "G")
  fragmentos_resultantes <- list(secuencia_adn)
  
  for (i in seq_along(nucleotidos)) {
    N <- nucleotidos[i]
    hinfi <- glue("GA{N}TC")
    nuevos_fragmentos <- list()
    
    # Aplicar la enzima a cada fragmento en la lista actual
    for (fragmento in fragmentos_resultantes) {
      fragmentos_cortados <- simular_enzima_restriccion(fragmento, hinfi, corte_pos)
      nuevos_fragmentos <- c(nuevos_fragmentos, fragmentos_cortados)
    }
    
    # Actualizar la lista de fragmentos resultantes
    fragmentos_resultantes <- nuevos_fragmentos
  }
  
  return(fragmentos_resultantes)
}

#Dibujar las bandas del gel electroforetico de Agarosa del ADN fragmentado
electroforesis_gel_agarosa <- function(fragmentos_DNA){
  nro_fragmentos <- nchar(fragmentos_DNA)
  ry <- range(nro_fragmentos) * c(0.5, 1.5)
  par (bg = "black")
  plot(1, type="n", ylim = ry, log = "y", xlim = c(0,2))
  
  
  for(i in nro_fragmentos){
    segments(0.5, i, 1.5, i, lwd=5, col = "coral")
  }
  
  axis(2, at = nro_fragmentos, labels = nro_fragmentos, pos = 0.4, las = 1, col = "white", col.axis = "white")
}

# Función para leer la secuencia desde un archivo .fasta
obtener_secuencia_en_fasta <- function(ruta_fasta) {
  secuencias <- readDNAStringSet(ruta_fasta)
  secuencia_adn <- as.character(secuencias[1]) # Suponiendo que solo hay una secuencia en el archivo
  return(secuencia_adn)
}

# Ejemplo de uso N°1 (Sin archivo .fasta - Secuencia directa)
#secuencia_adn <- "AAAAAAAGATTCAAAAGATTCAAAAAGACTCAAAAAAGAGTCGGGGGGCCCCCCAAAAAGAATCGGGGGCCCCCCCCCC"
#corte_pos <- 3  # La enzima corta después del tercer nucleótido en el sitio de restricción
#fragmentos_resultantes <- simular_corte_hinfi(secuencia_adn, corte_pos)



# Ejemplo de uso N°2 (con Archivo .fasta)
ruta_fasta <- "Ruta del archivo .fasta"
secuencia_adn <- obtener_secuencia_en_fasta(ruta_fasta)
corte_pos <- 3

fragmentos_resultantes <- simular_corte_hinfi(secuencia_adn, corte_pos)

# Imprimir los fragmentos resultantes
print(fragmentos_resultantes)
# Dibujar los fragmentos como bandas en el gel
electroforesis_gel_agarosa(fragmentos_resultantes)
