################################################################################

shuffleAndExtract <- function(inputFastaFile,
                              numberOfoutputSeqs,
                              lengthOutputSeqs,
                              initialPos=1,
                              outputFileName = paste(inputFastaFile,"output",
                                                     "fasta",sep=".")
                              ){

    ################################################################################
    #    handleFastaDatasets.R   
    #
    #    This function in R is designed to open a fasta file dataset, shuffle the 
    #    sequences and extract the desired sequences wanted by the user to generate 
    #    a new dataset of fixed size (number of required sequences) and with the 
    #    same length for each sequence. 
    #
    #    Author: Benjamin Tovar
    #    Post: http://btovar.com/2015/03/handling-fasta/
    #    Date: 22/JUNE/2012
    #
    ################################################################################

    #    run example:
    #    source("handleFastaDatasets.R")
    #    shuffleAndExtract("example.fasta",1000,200)

    # inputFastaFile: name of the input fasta file
    # numberOfoutputSeqs: number of desired sequences in the output file 
    # lengthOutputSeqs: fixed length of every sequence in the output file
    # initialPos: Position where the new window sizing will begin, default = 1
    # outputFileName: name of the output file, by default will be (e.g):
    #                 "inputFastaFile.fasta.output.fasta"

    cat("*** Starting computations |",date(),"****\n")    
        
    # Load the seqinr package
    require(seqinr)

    # Load the large seq and not shuffled dataset
    inputSeqs <- read.fasta(inputFastaFile)

    cat("\tProcessing for",length(inputSeqs),"sequences |",date(),"\n")    

    # Extract the length of every sequence and store the results in a vector.
    inputSeqsSize <- rep(NA,length(inputSeqs))
    for(i in 1:length(inputSeqs)){ 
        inputSeqsSize[i] <- summary(inputSeqs[[i]])$length
    }

    # Extract the index of the sequences which are longer than threshold
    inputSeqsLongSeqIndex <- which(inputSeqsSize > (lengthOutputSeqs+initialPos))

    # randomly pick numberOfoutputSeqs indexes that will be used to create 
    # the output dataset 
    inputSeqsIndex <- sample(inputSeqsLongSeqIndex,numberOfoutputSeqs,rep=F)
    
    # Store the Fasta header of each selected sequence in a vector
    inputSeqsIndexNames <- rep(NA,numberOfoutputSeqs)
    
    # create output object 
    outputSeqs <- list()
    for(i in 1:numberOfoutputSeqs){
        # Extract the fasta headers 
        inputSeqsIndexNames[i] <- attr(inputSeqs[[inputSeqsIndex[i]]],"name")
        # Extract the sequence
        outputSeqs[[i]] <- inputSeqs[[inputSeqsIndex[i]]][initialPos:((initialPos+lengthOutputSeqs)-1)]
    }
 
    # Export the sequences in a new fasta file 
    write.fasta(outputSeqs,inputSeqsIndexNames,outputFileName)
    
    cat("\n*** DONE |",date(),"****\n")
}


