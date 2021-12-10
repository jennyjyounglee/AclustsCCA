# sample annotation file creation code
library("IlluminaHumanMethylation450kanno.ilmn12.hg19")
annot <- minfi::getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
annot <- data.table(data.frame(annot))
annot <- annot[!chr%in%c("chrY","chrX"),]
annot <- annot[,list("IlmnID"=Name,"Coordinate_37"=pos,"CHR"=as.numeric(gsub("chr", "", chr)),Islands_Name,Relation_to_Island,UCSC_RefGene_Name)]
annot <- annot[CHR==7,]
annot <- annot[order(CHR,Coordinate_37),]

usethis::use_data(annot, compress = "xz", overwrite=T)
