setwd("your/working/dir")
#Filter out indels. They are marked as "1" in the list of SNPs.
library(dplyr)
list <- read.table(file="your_calls.file.txt",header=F)

list.snps <- list[which(list[,6]!=1),]

# Rename and reformat the data. 3 is chromosome number and 10934455 is the 
#position of the first nucleotide in studied short interval (Chilli Pepper).
list.snps$V1 <- 3
list.snps$V2 <- list.snps[,2] + 10934455
#If needed, filter out SNPs with lower quality (column 5).
#list.snps <- list.snps %>%
  #filter(V5 > 100)

all.split <- NULL
for(a in 1:length(list.snps[,1])){
  call <- as.character(list.snps[a,8])
  split <- strsplit(call,",")[[1]]
  all.split <- rbind(all.split,split)
}
all.split <- as.data.frame(all.split)
all.split <- sapply(all.split,as.numeric)

ref.count <-all.split[,1] + all.split[,2]
var.count <- all.split[,3] + all.split[,4]
list.snps.count <- cbind(list.snps[,1:7],ref.count,var.count)
blank <- rep(0,length(list.snps.count[,1]))
list.snps.count <- cbind(3, list.snps.count[,2:3], blank, list.snps.count[,4], blank)
write.table(list.snps.count, file="your_file.txt")

lib.nums <- seq(1,5)

for(a in 1:length(lib.nums)){
  
  calls <- read.table(file=paste("your_indiv_file", lib.nums[a], ".targcalls.txt",sep=""))
  
  calls$V1 <- "Chr3"
  calls$V2 <- calls$V2 + 10934455  
  
  calls <- calls[which(calls[,6] != 1),]
  
  values <- matrix(calls[,8], ncol=1)
  split <- apply(values, 1, function(x) strsplit(x, split=","))
  split.values <- matrix(unlist(split),ncol=4,byrow=T)
  split.values <- as.data.frame(split.values)
  split.values <- sapply(split.values, as.numeric)
  ref.count <- split.values[,1] + split.values[,2]
  var.count <- split.values[,3] + split.values[,4]
  calls.count <- cbind(calls[,1:7],ref.count,var.count)
  
  list.snps.count <- read.table("your_file.txt") 
  chp.complete <- NULL
  
  chp.coords <- list.snps.count[,2]
  calls.coords <- calls.count[,2]
  chp.match <- chp.coords %in% calls.coords
  calls.match <- calls.coords %in% chp.coords
  
  list.snps.count[which(chp.match==T),4] <- calls.count[which(calls.match==T),8]
  list.snps.count[which(chp.match==T),6] <- calls.count[which(calls.match==T),9]
  chp.complete <- rbind(chp.complete,list.snps.count)
  
  write.table(chp.complete,file=paste("your_indiv_file.", lib.nums[a] ,".allele_targ.complete.txt",sep=""),
              col.names=F,row.names=F,quote=F,sep="\t")
  print(paste("Sample ",as.character(a)," complete!", sep = ""))
}

#Load libraries
library(zoo)
library(ggplot2)
library(gridExtra)
library(lattice)
library(grid)
library(gtable)
library(xlsx)
library(ggrepel)
#Set coordinates of genes within interval and primer binding sites (optional)

HB34pos <- c(10940598,1,10941536,5)
lncRNApos <- c(10951458, 1, 10951757, 5)
amplicon1 <- c(10934673, 10943436)
amplicon2 <- c(10941447, 10951659)
amplicon3 <- c(10950600, 10960743)
samples <-  seq(1, 5)
#Start pdf and plotting genotypes
pdf("ESIL.pdf", onefile = T, height = 8, width = 16)

for(k in 1:length(samples)){
  file <- read.table(file=paste("your_indiv_file.",samples[k],".allele_targ.complete.txt", sep = ""))
  file$total_reads <- file$V4 + file$V6
  #Filter out SNPs with reads less than 30. 
  file <- file[file$total_reads >30,]
  if(nrow(file) > 0) {
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
    roll_vector
    
    call_SNP <- NULL
    max(adjLER1$adjLER1) - roll_vector
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
        call_SNP[i] <- "C24"
      } else if((max(adjLER3$adjLER3 - roll_vector[i])) >= 20 & !is.na(roll_vector[i]) & (max(adjLER3$adjLER3) - roll_vector[i]) < 75 & !is.na(roll_vector[i])){
        call_SNP[i] <- "HET"
      } else if((max(adjLER3$adjLER3 - roll_vector[i])) > 75 & !is.na(roll_vector[i])){
        call_SNP[i] <- "COL"
      }
      else if(is.na(roll_vector[i]) == T){
        next
      }}
    call_SNP
    
    if(length(call_SNP) == length(roll_vector)){
      genotype <- call_SNP
      Positions <- c(adjLER1$V2,adjLER2$V2,adjLER3$V2)
      ReadsLER <- c(vectorLER)
      ReadsCOL <- c(vectorCOL)
      snpLER <- c(adjLER1$V5,adjLER2$V5,adjLER3$V5)
      snpCOL <- c(adjCOL1$V3,adjCOL2$V3,adjCOL3$V3)
      logreads <- c(file$logreads)
      CO <-  data.frame(Positions,ReadsLER,ReadsCOL,snpLER,snpCOL,genotype)
    } else if(length(call_SNP) != length(roll_vector)){
      next
    }
    #Write CO file with data
    write.xlsx(x = CO, file = paste("your_indiv_file", samples[k], ".crossover.xlsx", sep=""),sheetName = as.character(samples[k]) , row.names = T, col.names = T, append = F) 
    ChP_length <- seq(10934455, 10934455 + 26430, 26431 / nrow(CO))
    #plot
    s <- ggplot(data = CO, aes(x = Positions)) +  
      geom_line(aes(y = vectorLER), col = "red", size = 0.8) + geom_point(aes(y = vectorLER), col = "red", size = 1.5) + 
      geom_line(aes(y = vectorCOL), col = "blue", size = 0.8) + geom_point(aes(y = vectorCOL), col = "blue", size = 1.5) + 
      #geom_segment(aes(x = amplicon1[1], y = 0, xend = amplicon1[1], yend = Inf), linetype = "dashed", col = "chocolate2", size = 0.35) +
      #geom_segment(aes(x = amplicon1[2], y = 0, xend = amplicon1[2], yend = Inf), linetype = "dashed", col = "chocolate2", size = 0.35) +
      #geom_segment(aes(x = amplicon2[1], y = 0, xend = amplicon2[1], yend = Inf), linetype = "dashed", col = "green", size = 0.35) +
      #geom_segment(aes(x = amplicon2[2], y = 0, xend = amplicon2[2], yend = Inf), linetype = "dashed", col = "green", size = 0.35) +
      #geom_segment(aes(x = amplicon3[1], y = 0, xend = amplicon3[1], yend = Inf), linetype = "dashed", col = "magenta", size = 0.35) +
      #geom_segment(aes(x = amplicon3[2], y = 0, xend = amplicon3[2], yend = Inf), linetype = "dashed", col = "magenta", size = 0.35) +
      theme_minimal() +
      #geom_rect(aes(xmin = HB34pos[1], xmax = HB34pos[3], ymin = HB34pos[2], ymax = HB34pos[4]), fill = "#F9E79F", color = "NA") +
      #geom_rect(aes(xmin = lncRNApos[1], xmax = lncRNApos[3], ymin = lncRNApos[2], ymax = lncRNApos[4]), fill = "#F9E79F", color = "NA") + 
      theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                         panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) 
    s <- s + ggtitle(paste(as.character(samples[k]), "Number of reads:", as.character(total_reads))) + labs(y = "Percent of reads", x = "ChP SNPs") + 
      theme(plot.title = element_text(hjust = 0.5)) 
    s <- s + coord_cartesian(xlim = c(min(ChP_length), max(ChP_length)),ylim = c(5,100))   
    #s <- s + ggrepel::geom_text_repel(data = CO, y = ReadsLER, label = as.character(Positions), size = 1, col = "purple", segment.size = 0.5, segment.alpha = 0.5, segment.color = "purple")
    
    plot(s)
  }
  else next
}
dev.off()

