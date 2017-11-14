# Functions for the RADTools package


#
# 1. Read and process FASTA File
#
process_fasta <- function(fasta_path){
  # Function to read FASTA and extract sequences. 
  # It converts multi-line sequences into single-line sequences.
  # Will output vector object containing all sequences.
  # Function by Angel G. Rivera-Colon
  
  fa <- file(fasta_path, 'r')
  all_seq <- c()      # Final vector containing all sequences in the reference
  loc_seq <- c()      # Local vector containing per-chromosome/scaffold sequences
  while( TRUE ){
    line <- readLines(fa, n=1)
    # break at EOF and output last sequence set
    if( length(line) == 0 ){
      # Process last sequence set
      all_seq <- c(all_seq, paste(loc_seq, collapse=''))
      break
    }
    # if line is a FASTA header (first character is ">")
    if(strsplit(line, split='')[[1]][1] == '>'){
      # And if not the first header read
      if(length(loc_seq) != 0){
        # process previous sequence set
        # append all sequence lines into a single long sequence
        all_seq <- c(all_seq, paste(loc_seq, collapse=''))
        # clear for next sequence set
        loc_seq <- c()
      }
    } else {
      # if it is a sequence line just append to the local sequence vector
      loc_seq <- c(loc_seq, line)
    }
  }
  close(fa)    # close conection 
  return(all_seq)  # output final sequence vector
}


#
# 2. Find cutsites in reference sequence
#
find_cuts <- function(sequences, cutsite){
  # Function to find restriction cutsites in a sequence of interest
  # For restriction site 'cutsite', find all the position is appears on sequence 'sequence'. 
  # Works on output of `process_fasta` function, vector containing multiple sequences.  
  # Returns a vector containing all the positions (indexes) of the cutsites in the sequence of interest. 
  # Final output is a list containing the vector positions for each of the sequences of interest. 
  # Function by Angel G. Rivera-Colon
  
  final_cut <- list()           # final output list, contains all cut vectors
  idx <- 1                      # index of final_cut list
  for(seq in sequences){
    # Find string `cutsite` in string `seq`
    pos <- gregexpr(cutsite, seq, fixed=T, useBytes=T)
    # Extract the position vector from regex object
    cuts <- pos[[1]][1:length(pos[[1]])]
    # Append it to the final list
    final_cut[[idx]] <- cuts
    idx <- idx + 1
  }
  return(final_cut)
}



#
# 3. Function to calculate inter-cutsite distance distribution
#
cutsite_distance <- function(cut_list){
  # Function to calculate inter-cutsite distance (distance from one cutsite to the other) distribution.
  # Input file is a list of vector containing the per chromosome/scaffold cutsite positions (output of `find_cuts`).
  # Output is a single vector containing the inter-cutsite distances for the whole genome. 
  # Output is compatible with other R tools like `summary()`, `hist()`, etc. 
  # Function by Angel G. Rivera-Colon
  
  distances_all <- c()                # Vector containing all calculated distances
  # Loop throught the chromosomes...
  for(cuts in cut_list){
    # Position of previous cutsite, defaults to 0 at the start of each chromosome. 
    prev <- 0
    # Loop through the cuts in each chromosome...
    for(cut in cuts){
      distance <- cut - prev
      distances_all <- c(distances_all, distance)
      prev <- cut
    }
  }
  return(distances_all)
}







# Define sequencing machines
hiseq2500 <- c(2.2e8, 3.1e8, 4e8)
names(hiseq2500) <- c('Low', 'Med', 'Hi')

hiseq4000 <- c(5.0e8, 7.5e8, 1.0e9)
names(hiseq4000) <- c('Low', 'Med', 'Hi')


Per_Sample_Coverage <- function(num_cutsites,                  #Number of cutsites present in your genome
                                num_samples,                   #Number of total samples you want to sequence
                                sequencing_machine = NULL){    #Illumina sequencing platform you want to use
  #This function calculates the predicted per sample
  #coverage from RADseq given the number of cutsites
  #in the genome, number of samples, and the type
  #of illumina machine
  #Function by Kira Long
  if(is.null(sequencing_machine)){                    #If no machine is provided, function will default hiseq2500
    sequencing_machine <- "hiseq2500"
  }
  if(sequencing_machine == "hiseq2500"){              #This if statement determines the number of reads
    ReadsPerLane <- hiseq2500
	#per lane in the desired sequencer with the range
  }else if(sequencing_machine == "hiseq4000"){        #of low, medium, or high total reads per lane
    ReadsPerLane <- hiseq4000
  }
  RADtags <- num_cutsites*2                           #Determines your number of RADtags
  DNAseqs <- num_samples * RADtags                    #Determines how many DNA sequences you will get from all samples
  coverage <- ReadsPerLane/DNAseqs                    #Calculates your coverage
  return(print(coverage, digits = 4))
}


Samples_Per_Lane <- function(num_cutsites, desired_coverage = NULL, sequencing_machine = NULL){
  #This function calculates the estimate number of samples you can put in 1 illumina lane to
  #get a desired, set coverage based on the number of cutsites in your genome and the illumina 
  #machine used to sequence
  #Function by Kira Long
  if(is.null(desired_coverage)){                      #If no desired coverage is provided, function will default to 30
    desired_coverage <- 30                            #as 30 is minimum to identify heterozygotes
  }
  if(is.null(sequencing_machine)){                    #If no machine is provided, function will default hiseq2500
    sequencing_machine <- "hiseq2500"
  }
  if(sequencing_machine == "hiseq2500"){              #This if statement determines the number of reads
    ReadsPerLane <- hiseq2500                         #per lane in the desired sequencer with the range
  }else if(sequencing_machine == "hiseq4000"){        #of low, medium, or high total reads per lane
    ReadsPerLane <- hiseq4000
  }
  RADtags <- num_cutsites*2
  min_DNA_seqs <- ReadsPerLane/desired_coverage
  estimated_num_samples <- min_DNA_seqs/RADtags
  return(round(estimated_num_samples))
}
