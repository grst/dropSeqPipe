library(stringdist)
reads_per_barcode = read.csv(snakemake@input$reads_per_barcode, sep = '\t')


reference = read.csv(snakemake@input$reference, header = FALSE, stringsAsFactors = FALSE)$V1

reads_per_barcode_filtered = subset(reads_per_barcode, !reads_per_barcode$TAG.XC %in% reference)

total = sum(reads_per_barcode_filtered[,1])

y_raw=cumsum(reads_per_barcode_filtered[,1])
y=(y_raw/total)

#90 percent of total reads are in X numbr of barcodes
id = sum(y < 0.9)

cbc_dist = stringdistmatrix(reference,reads_per_barcode_filtered[1:id,2], useNames = TRUE)
kept_barcodes = c()
max_dist=1
for (dist in 1:max_dist){
  ids = which(apply(cbc_dist, 2, function(x) sum(x==1)) == 1)
  kept_barcodes = c(kept_barcodes, colnames(cbc_dist)[ids])
  cbc_dist = cbc_dist[,-ids]
}
if (identical(kept_barcodes, character(0))){
  write.table(file=snakemake@output[[1]], data.frame(neighbour='None', reference='None'), row.names = FALSE, sep = '\t', quote = FALSE)
  }else{

new = stringdistmatrix(reference, kept_barcodes, useNames = TRUE)
mapping = data.frame(matrix(ncol=2, nrow=ncol(new)))
colnames(mapping) = c('neighbour','reference')
start=1
for (i in 1:length(reference)){
  neighbours = names(which(new[i,] == 1))
  if(length(neighbours) == 0){next}
  end = start + length(neighbours)-1
  temp = data.frame(neighbour=neighbours, reference=cbind(rep(reference[i],length(neighbours))), stringsAsFactors = FALSE)
  mapping[start:end,] = temp
  start=start+length(neighbours)
}

write.table(mapping, snakemake@output[[1]], row.names = FALSE, sep = '\t', quote = FALSE)
}