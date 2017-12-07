# Functions for the RADTools package
# BitBucket repository: https://bitbucket.org/angelgr2/radseq_tools


#
# 1. Read and process FASTA File
#
process_fasta <- function(fasta_path, 
                          min_seq_len=NULL){
  # Function to read FASTA and extract sequences. 
  # It converts multi-line sequences into single-line sequences.
  # Will output vector object containing all sequences.
  # Function by Angel G. Rivera-Colon
  
  # If no minimum sequence length is specified, use default value (5000)
  if(is.null(min_seq_len)){                    
    min_seq_len <- 5000
  }
  
  fa <- file(fasta_path, 'r')
  all_seq <- c()      # Final vector containing all sequences in the reference
  loc_seq <- c()      # Local vector containing per-chromosome/scaffold sequences
  long_seq <- NULL    # Variable containing the merged sequences
  while( TRUE ){
    line <- readLines(fa, n=1)
    # break at EOF and output last sequence set
    if( length(line) == 0 ){
      # Process last sequence set
      
      # Append all sequence lines into a single long sequence
      long_seq <- paste(loc_seq, collapse='')
      
      # If sequence is long enough, add to `all_seq` vector
      if(nchar(long_seq) >= min_seq_len){
        all_seq <- c(all_seq, long_seq)
      }
      break
      
    }
    # if line is a FASTA header (first character is ">")
    if(strsplit(line, split='')[[1]][1] == '>'){
      # And if not the first header read
      if(length(loc_seq) != 0){
        # Process previous sequence set
        
        # Append all sequence lines into a single long sequence
        long_seq <- paste(loc_seq, collapse='')
        
        # If sequence is long enough, add to `all_seq` vector
        if(nchar(long_seq) >= min_seq_len){
          all_seq <- c(all_seq, long_seq)
        }
        
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

# Vector containing commonly used restriction enzymes
renz <- c(
  pstI = 'CTGCAG',
  sbfI = 'CCTGCAGG',
  ecoRI = 'GAATTC',
  hindIII = 'AAGCTT'
)

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
    # If cut object is not `-1` (no objects found), append it to the final list
    if(cuts[1] != -1){
      final_cut[[idx]] <- cuts
    } else {
      final_cut[[idx]] <- NA
    }
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
    # Add an initial zero to the `cuts` vector, so distance of first cutsite can be calculated. 
    cuts <- c(0, cuts)
    # Calculate distance between cut sites using `diff()`. Omit missing values (for sequences without cuts)
    distances <- diff(na.omit(cuts))
    # Add `distances` vector (this chromosome) to `distances_all` (all chromosomes)
    distances_all <- c(distances_all, distances)
  }
  return(distances_all)
}

#
#4. Find number of cutsites 
#
number_cutsites <- function(cutsite_list){ #Cutsite_list is an object product of the find_cuts() function
#This function calculates the total number of cutsites in the genome
#Function by Kira Long
  count <- 0
  for (i in 1:length(cutsite_list)){
    #Skip cut vector if it contains NA i.e. a cutsite was not found in the sequence
    if(is.na(cutsite_list[[i]][1]) == FALSE){
      count <- count + length(cutsite_list[[i]])
    }
  }
  #Returns a single numeric value, which is the sum of all cutsites
  return(count)
}

#Defining sequencing machines for functions 5 and 6
hiseq2500 <- c(2.2e8, 3.1e8, 4e8)
names(hiseq2500) <- c('Low', 'Med', 'Hi')

hiseq4000 <- c(5.0e8, 7.5e8, 1.0e9)
names(hiseq4000) <- c('Low', 'Med', 'Hi')

#
#5. Find the coverage per sample
#
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


#
#6. Find the number of samples that can fit in one illumina lane
#
Samples_Per_Lane <- function(num_cutsites,
			     desired_coverage = NULL,
			     sequencing_machine = NULL){
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


#
# 7. Find the number of DNA sequence reads
#
DNA_Reads_Per_Lane <- function(num_cutsites, 
                               num_samples,
                               desired_coverage = NULL){
  #This function calculates the estimated number of DNA reads that will be generated based
  #on the number of cutsites, number of samples, and desired coverage.
  #Function by Kira Long
  if(is.null(desired_coverage)){    #If no desired coverage is provided, function will default to 30
    desired_coverage <- 30          #as 30 is minimum to identify heterozygotes
  }
  RADtags <- num_cutsites*2         #Determines your number of RADtags
  #Determine how many DNA sequences you will get from all samples
  DNA_reads <- num_samples * RADtags * desired_coverage
  return(print(DNA_reads, digits=3)) #Returns number of reads
}
