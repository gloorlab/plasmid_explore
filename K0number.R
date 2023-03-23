input <- c('1-T1', '2-T4')


get_kegg <- function(sample, verbose=FALSE){
      dir <- paste('~/Documents/thesis/data/', sample,'/', sep='')

      x <- list.files(dir, pattern='[0-9].tsv') # tsv file only
        # move up.files to a simple variable for convenience
        
        annot <- matrix(data=NA, nrow=length(x), ncol=5)
        len <- vector()
        keggFound <- vector()
        ORFS <- vector()
        
        for(i in 1:length(x)){
            if(verbose) print(paste(dir,x[i], sep=''))
            gb <- read.table(paste(dir,x[i], sep=''), comment.char="#", sep='\t', quote='')
            annot[i,1] <- dim(gb)[1] # total cds
            annot[i,2] <- length(grep('PFAM', gb[,9])) # PFAM annotated ORFS in column 9
            kolines <- grep('KEGG', gb[,9])
            for(k in kolines){ keggFound <- c(keggFound, gsub(".*KEGG:(K\\d+).+", "\\1", gb[k,9]) )}
            annot[i,4] <- length(kolines) # KEGG annotated ORFS in column 9
            annot[i,3] <- annot[i,2]/annot[i,1] # proportion n.annot / n.rows
            annot[i,5] <- annot[i,4]/annot[i,1] # proportion n.annot / n.rows
            len <- c(len, gb[,4] - gb[,3])
        }
        output <- list('keggFound'=unique(keggFound), 'nPlasmids' = length(x), 'nORFs'=sum(annot[,1]))
	return(output)
}
unique(keggFound) # unique KO found
x # number of plasmids
sum(annot[,1]) # number of ORFS