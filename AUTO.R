library(dplyr)
library(zoo)
library(ggplot2)
library(gridExtra)
library(lattice)
library(grid)
library(gtable)
library(xlsx)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  print("List of .allele_targ.complete.txt files is necessary!")
  stop("Requires command line argument.")}


list1 <- read.table(args[1], stringsAsFactors = FALSE)
HB34pos <- c(10940598,1,10941536,5)
lncRNApos <- c(10951458, 1, 10951757, 5)
amplicon1 <- c(10934673, 10943436)
amplicon2 <- c(10941447, 10951659)
amplicon3 <- c(10950600, 10960743)

pdf(paste(strsplit(list1[1,1], "\\.")[[1]][1],".pdf",sep=""), onefile = T, height = 8, width = 16)
for(k in 1:nrow(list1)){
  file_name <- list1[k,1]
  if(file.exists(file_name)){
    if(file.info(file_name)$size > 0){
    file <- read.table(file_name)
      file$total_reads <- file$V4 + file$V6
      #Filter out SNPs with reads less than 30. 
      file <- file[file$total_reads >100,]
      #Calculate percentage of reads 
      file$percV4 <- file$V4 * 100 / file$total_reads
      file$percV6 <- file$V6 * 100 / file$total_reads
      total_reads <- sum(file$total_reads)
      #Calculate logarithm of total number of reads
      file$logreads <- log(file$total_reads)
      #Subset interval into amplicons based on position of primers. You can adjust
      #the percentage if you observe some amplification bias in your long-range PCR.
      adjLER1 <- subset(file, file[,2] < amplicon1[2])
      adjLER1$adjLER1 <- adjLER1$percV6 * 1
      adjLER2 <- subset(file, (file[,2] > amplicon1[2]) & (file[,2] < amplicon2[2]))
      adjLER2$adjLER2 <- adjLER2$percV6 * 1
      adjLER3 <- subset(file, (file[,2] > amplicon2[2]) & (file[,2]< amplicon3[2]))
      adjLER3$adjLER3 <- adjLER3$percV6 * 1
      
      adjCOL1 <- subset(file, file[,2] < amplicon1[2])
      adjCOL1$adjCOL1 <- adjCOL1$percV4 * 1
      adjCOL2 <- subset(file, (file[,2] > amplicon1[2]) & (file[,2] < amplicon2[2]))
      adjCOL2$adjCOL2 <- adjCOL2$percV4 * 1
      adjCOL3 <- subset(file, (file[,2] > amplicon2[2]) & (file[,2]< amplicon3[2]))
      adjCOL3$adjCOL3 <- adjCOL3$percV4 * 1
      
      vectorLER <- c(adjLER1$adjLER1,adjLER2$adjLER2,adjLER3$adjLER3)
      vectorCOL <- c(adjCOL1$adjCOL1,adjCOL2$adjCOL2,adjCOL3$adjCOL3)
      roll_vector <- rollmean(vectorLER, k = 3, fill = mean(vectorLER[1:3]))
      call_SNP <- NULL
      if(length(roll_vector) != 0){
        for(i in 1:length(roll_vector)){
          if((max(adjLER1$adjLER1 - roll_vector[i])) <= 20 & !is.na(roll_vector[i])){
            call_SNP[i] <- "LER"
          } else if((max(adjLER1$adjLER1 - roll_vector[i])) >= 20 & !is.na(roll_vector[i]) & (max(adjLER1$adjLER1) - roll_vector[i]) < 75 & !is.na(roll_vector[i])){
            call_SNP[i] <- "HET"
          } else if((max(adjLER1$adjLER1 - roll_vector[i])) > 75 & !is.na(roll_vector[i])){
            call_SNP[i] <- "COL"
          } else if((max(adjLER2$adjLER2 - roll_vector[i])) <= 20 & !is.na(roll_vector[i])){
            call_SNP[i] <- "LER"
          } else if((max(adjLER2$adjLER2 - roll_vector[i])) >= 20 & !is.na(roll_vector[i]) & (max(adjLER2$adjLER2) - roll_vector[i]) < 75 & !is.na(roll_vector[i])){
            call_SNP[i] <- "HET"
          } else if((max(adjLER2$adjLER2 - roll_vector[i])) > 75 & !is.na(roll_vector[i])){
            call_SNP[i] <- "LER"
          } else if((max(adjLER3$adjLER3 - roll_vector[i])) <= 20 & !is.na(roll_vector[i])){
            call_SNP[i] <- "COL"
          } else if((max(adjLER3$adjLER3 - roll_vector[i])) >= 20 & !is.na(roll_vector[i]) & (max(adjLER3$adjLER3) - roll_vector[i]) < 75 & !is.na(roll_vector[i])){
            call_SNP[i] <- "HET"
          } else if((max(adjLER3$adjLER3 - roll_vector[i])) > 75 & !is.na(roll_vector[i])){
            call_SNP[i] <- "COL"
          }
          else if(is.na(roll_vector[i]) == T){
            next
          }}
      
      if(length(call_SNP) == length(roll_vector)){
        genotype <- call_SNP
        Positions <- c(adjLER1$V2,adjLER2$V2,adjLER3$V2)
        ReadsLER <- c(vectorLER)
        ReadsCOL <- c(vectorCOL)
        snpLER <- c(adjLER1$V5,adjLER2$V5,adjLER3$V5)
        snpCOL <- c(adjCOL1$V3,adjCOL2$V3,adjCOL3$V3)
        logreads <- c(file$logreads)
        CO <-  data.frame(Positions,ReadsLER,ReadsCOL,snpLER,snpCOL,genotype)
      } 
      else if(length(call_SNP) != length(roll_vector)){
      next
      }
      parts <- strsplit(file_name, "\\.")[[1]]
      name <- parts[c(1, 2)]
      name_combined <- paste(name, collapse = ".")
      write.xlsx(x = CO, file = paste(name_combined, ".crossover.xlsx", sep=""), 
                  row.names = T, col.names = T, append = F) 
      ChP_length <- seq(10934455, 10934455 + 26430, 26431 / nrow(CO))
      #plot
      s <- ggplot(data = CO, aes(x = Positions)) +  
        geom_line(aes(y = vectorLER), col = "red", size = 0.8) + geom_point(aes(y = vectorLER), col = "red", size = 1.5) + 
        geom_line(aes(y = vectorCOL), col = "blue", size = 0.8) + geom_point(aes(y = vectorCOL), col = "blue", size = 1.5) + 
        theme_minimal() +
        geom_rect(aes(xmin = HB34pos[1], xmax = HB34pos[3], ymin = HB34pos[2], ymax = HB34pos[4]), fill = "#F9E79F", color = "NA") +
        geom_rect(aes(xmin = lncRNApos[1], xmax = lncRNApos[3], ymin = lncRNApos[2], ymax = lncRNApos[4]), fill = "#F9E79F", color = "NA") + 
        theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) 
      s <- s + ggtitle(paste(as.character(name_combined), "Number of reads:", as.character(total_reads))) + labs(y = "Percent of reads", x = "ChP SNPs") + 
        theme(plot.title = element_text(hjust = 0.5)) 
      s <- s + coord_cartesian(xlim = c(min(ChP_length), max(ChP_length)),ylim = c(5,100))   
      
      plot(s)
      }else next
    } else next
  } else next
}  
dev.off()

