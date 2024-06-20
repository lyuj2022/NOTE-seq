library(dplyr)

args <- commandArgs(T)
samplename <- args[1]
cutoff <- 0.6

blast <- read.table(paste0(samplename,'_both_strand_available_site.tsv'), fill=T)[,c(1,2,6,8,9)]
info <- read.table(paste0(samplename,'_filtered_barcode_and_gene.txt'), row.names = 1)##### name + '\t' + cb + umi + '\t' + gene + '\n'
outname <- paste0(samplename,'_count.txt')

blast <- blast[which(blast$V1 %in% rownames(info)),]
blast_add_info <- cbind(blast,info[blast$V1,])
rm(blast)
colnames(blast_add_info) <- c('Name','Flag','READ BASE','REF POS','REF BASE','label','Gene')


blast_group_by_umi <- split(blast_add_info,blast_add_info$label)
mutation_count <- lapply(blast_group_by_umi, function(x){
  flag <- x[1,2]
  able_site <- length(unique(x$`REF POS`))
  freq_count <- table(x$`READ BASE`,x$`REF POS`)
  freq_count_mtx <- as.matrix(t(freq_count))
  freq <- prop.table(freq_count_mtx,1)
  freq <- as.data.frame.array(freq)
  base <- colnames(freq)
  if(flag == 0){
    if('C' %in% base){
      freq_site_count <- freq[which(freq$C >= cutoff),]
      mut_num <- dim(freq_site_count)[1]
    }
    else{
    mut_num <- 0
    }
  }
  else{
    if('G' %in% base){
      freq_site_count <- freq[which(freq$G >= cutoff),]
      mut_num <- dim(freq_site_count)[1]
    }
    else{
    mut_num <- 0
    }
  }
  return(c(able_site, mut_num))
})

mutation_count <- do.call(rbind,mutation_count)
mutation_count <- as.data.frame(mutation_count)

rownames(info) <- info$V2
gene <- info[rownames(mutation_count),]
label <- rownames(mutation_count)
barcode <- substr(label, 1, 18)
mutation_count$gene <- gene
mutation_count$cell <- barcode

write.table(mutation_count,outname,quote = F,col.names = F)

