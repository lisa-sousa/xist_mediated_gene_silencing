
####input files
yang_2010_file = "/project/lncrna/Xist/data/annotation_files/escapees/metadata/escapees_2010_yang.txt"
splinter_2011_file = "/project/lncrna/Xist/data/annotation_files/escapees/metadata/escapees_2011_splinter.txt"
calabrese_2012_file = "/project/lncrna/Xist/data/annotation_files/escapees/metadata/escapees_2012_calabrese.txt"
wu_2014_file = "/project/lncrna/Xist/data/annotation_files/escapees/metadata/escapees_2014_wu.txt"
berlecht_2015_file = "/project/lncrna/Xist/data/annotation_files/escapees/metadata/escapees_2015_berletch.txt"
marks_2015_file = "/project/lncrna/Xist/data/annotation_files/escapees/metadata/escapees_2015_marks.txt"
andergassen_2017_random_file = "/project/lncrna/Xist/data/annotation_files/escapees/metadata/escapees_2017_andergassen_random.txt"
andergassen_2017_imprinted_file = "/project/lncrna/Xist/data/annotation_files/escapees/metadata/escapees_2017_andergassen_imprinted.txt"

####load escapee lists
yang_2010 = read.table(yang_2010_file)
splinter_2011 = read.table(splinter_2011_file)
calabrese_2012 = read.table(calabrese_2012_file)
wu_2014 = read.table(wu_2014_file)
berlecht_2015 = read.table(berlecht_2015_file)
marks_2015 = read.table(marks_2015_file)
andergassen_2017_random = read.table(andergassen_2017_random_file)
andergassen_2017_imprinted = read.table(andergassen_2017_imprinted_file)

####merge tables
yang_2010$V1 = tolower(yang_2010$V1)
splinter_2011$V1 = tolower(splinter_2011$V1)
calabrese_2012$V1 = tolower(calabrese_2012$V1)
wu_2014$V1 = tolower(wu_2014$V1)
berlecht_2015$V1 = tolower(berlecht_2015$V1)
marks_2015$V1 = tolower(marks_2015$V1)
andergassen_2017_random$V1 = tolower(andergassen_2017_random$V1)
andergassen_2017_imprinted$V1 = tolower(andergassen_2017_imprinted$V1)

escapees = merge(yang_2010,splinter_2011,by="V1",all=T)
escapees = merge(escapees,calabrese_2012,by="V1",all=T)
escapees = merge(escapees,wu_2014,by="V1",all=T)
escapees = merge(escapees,berlecht_2015,by="V1",all=T)
escapees = merge(escapees,marks_2015,by="V1",all=T)
escapees = merge(escapees,andergassen_2017_random,by="V1",all=T)
escapees = merge(escapees,andergassen_2017_imprinted,by="V1",all=T)
colnames(escapees) = c("Genes","1","2","3","4","5","6","7","8")
escapees = escapees[order(escapees$Genes),]

write.table(escapees, file = "/project/lncrna/Xist/data/annotation_files/escapees/escapees_literature.txt", col.names = T, row.names = F, quote = F, sep = "\t")
