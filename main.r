# Librerías
install.packages("seqinr")
install.packages("tidyverse")
install.packages("gridExtra")

install.packages("BiocManager")
install.packages("ade4")
install.packages("ape")
install.packages("viridis")
install.packages("Rtools")
BiocManager::install("Biostrings")
BiocManager::install("DECIPHER")
library(ggplot2)
library(gridExtra)
library(ade4)
library(ape)
library(DECIPHER)
library(Biostrings)
library(viridis)
library(seqinr)
library(rtools)
if(!require('ade4')) {
  install.packages('ade4')
  library('ade4')
}

# Primero vamos a establecer los genomas
covid.variants <- c(
  "AY278489",
  "AY390556",
  "AY485277",
  "AY508724",
  "JX869059",
  "MN908947",
  "MN985325",
  "MT292571",
  "MT873893.1",
  "MW000351.1"
)

virus.sequences <- read.GenBank(covid.variants)

write.dna(virus.sequences, file="coronavirus_seqs.fasta", format="fasta")
virus.seq.not.align <- readDNAStringSet("coronavirus_seqs.fasta",
                                        format="fasta")


# 1
crear.grafica.longitud <- function(genomes) {
  # Creamos un dataframe vacio
  frame.size <- data.frame(matrix(NA, nrow = 10, ncol = 2))
  colnames(frame.size) <- c("Variantes", "Tamaño")
  
  # Ponemos las variantes
  frame.size["Variantes"] <- covid.variants
  
  # Calculamos los tamaños
  vector.size <- c()
  for(i in 1:length(genomes)) {
    fasta.sequence <- genomes[i]
    
    size <- nchar(fasta.sequence[[1]])
    vector.size <- c(vector.size, size)
  }
  
  # Ponemos los tamaños
  frame.size["Tamaño"] <- vector.size
  
  # Imprimimos el Dataframe
  print(frame.size)
  
  # Creamos la gráfica y la imprimimos
  plot.size.genomes <- ggplot(frame.size, aes_string(x="Variantes", y="Tamaño",
      fill="Variantes")) + geom_bar(stat="identity")
  print(plot.size.genomes)
}

crear.grafica.gc <- function(genomes) {
  # Creamos un dataframe vacio
  frame.size <- data.frame(matrix(NA, nrow = 10, ncol = 2))
  colnames(frame.size) <- c("Variantes", "GC")
  
  # Ponemos las variantes
  frame.size["Variantes"] <- covid.variants
  
  # Calculamos el porcentaje de GC
  vector.gc <- c()
  for(i in 1:length(genomes)) {
    fasta.sequence <- genomes[i]
    size <- nchar(fasta.sequence[[1]])
    size.gc <- 0
    
    for(j in 1:size) {
      character <- substr(fasta.sequence, j, j)
      
      if (character == "G" || character == "C") {
        size.gc <- size.gc + 1
      }
    }
    
    percentage <- size.gc / size * 100

    vector.gc <- c(vector.gc, percentage)
  }
  
  # Ponemos los tamaños
  frame.size["GC"] <- vector.gc
  
  # Imprimimos el Dataframe
  print(frame.size)
  
  # Creamos la gráfica y la imprimimos
  plot.gc.genomes <- ggplot(frame.size, aes_string(x="Variantes",
      y="GC", fill="Variantes")) + geom_bar(stat="identity")
  print(plot.gc.genomes)
}

crear.grafica.longitud(virus.seq.not.align)
crear.grafica.gc(virus.seq.not.align)



# 2
crear.grafica.agct <- function(genomes) {
  genomes.frame <- data.frame(
    "Nucleotides" = c("A", "G", "C", "T")
  )
  for(i in 1:length(genomes)) {
    fasta.sequence <- genomes[i]
    size <- nchar(fasta.sequence[[1]])
    nucleotides <- c(0, 0, 0, 0)
    
    for(j in 1:size) {
      character <- substr(fasta.sequence, j, j)
      
      if (character == "A") {
        nucleotides[1] <- nucleotides[1] + 1
      } else if (character == "G") {
        nucleotides[2] <- nucleotides[2] + 1
      } else if (character == "C") {
        nucleotides[3] <- nucleotides[3] + 1
      } else if (character == "T") {
        nucleotides[4] <- nucleotides[4] + 1
      }
    }
    
    # Data Frame
    composition.df =  c(nucleotides[1],
                        nucleotides[3],
                        nucleotides[2],
                        nucleotides[4]
    )
    
    genomes.frame[covid.variants[i]] <- unlist(composition.df)
  }
  
  print(genomes.frame)
  
  # Gráficas
  plot_count_virus_1 <- ggplot(genomes.frame, aes_string(x="Nucleotides",
    y=covid.variants[1], fill="Nucleotides")) + geom_bar(stat="identity")
  plot_count_virus_2 <- ggplot(genomes.frame, aes_string(x="Nucleotides",
    y=covid.variants[2], fill="Nucleotides")) + geom_bar(stat="identity")
  plot_count_virus_3 <- ggplot(genomes.frame, aes_string(x="Nucleotides",
    y=covid.variants[3], fill="Nucleotides")) + geom_bar(stat="identity")
  plot_count_virus_4 <- ggplot(genomes.frame, aes_string(x="Nucleotides",
    y=covid.variants[4], fill="Nucleotides")) + geom_bar(stat="identity")
  plot_count_virus_5 <- ggplot(genomes.frame, aes_string(x="Nucleotides",
    y=covid.variants[5], fill="Nucleotides")) + geom_bar(stat="identity")
  plot_count_virus_6 <- ggplot(genomes.frame, aes_string(x="Nucleotides",
    y=covid.variants[6], fill="Nucleotides")) + geom_bar(stat="identity")
  plot_count_virus_7 <- ggplot(genomes.frame, aes_string(x="Nucleotides",
    y=covid.variants[7], fill="Nucleotides")) + geom_bar(stat="identity")
  plot_count_virus_8 <- ggplot(genomes.frame, aes_string(x="Nucleotides", 
    y=covid.variants[8], fill="Nucleotides")) + geom_bar(stat="identity")
  plot_count_virus_9 <- ggplot(genomes.frame, aes_string(x="Nucleotides", 
    y=covid.variants[9], fill="Nucleotides")) + geom_bar(stat="identity")
  plot_count_virus_10 <- ggplot(genomes.frame, aes_string(x="Nucleotides", 
    y=covid.variants[10], fill="Nucleotides")) + geom_bar(stat="identity")
  
  plot.counts <- list(
    plot_count_virus_1,
    plot_count_virus_2,
    plot_count_virus_3,
    plot_count_virus_4,
    plot_count_virus_5,
    plot_count_virus_6,
    plot_count_virus_7,
    plot_count_virus_8,
    plot_count_virus_9,
    plot_count_virus_10
  )
  
  grid.arrange(grobs = plot.counts, ncol=2)
}

crear.grafica.agct(virus.seq.not.align)


# 6
# Orientar
virus.seq.not.align <- OrientNucleotides(virus.seq.not.align)

# Alinear
virus.seq.align <- AlignSeqs(virus.seq.not.align)

# Obtener el fasta
writeXStringSet(virus.seq.align, file="coronavirus_seq_align.fasta")
virus.aligned <- read.alignment("coronavirus_seq_align.fasta", format="fasta")

# Matriz
matriz.distancia <- dist.alignment(virus.aligned, matrix="similarity")
temp <- as.data.frame(as.matrix(matriz.distancia))

# Tabla de grises
table.paint(temp,clabel.row=.5, clabel.col=.5,cleg=0) + scale_color_viridis()

# 7
# Crear el árbol
virus.tree <- nj(matriz.distancia)

# Arreglarlo para su visualización
virus.tree <- ladderize(virus.tree)
plot(virus.tree)


