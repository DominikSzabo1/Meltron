suppressMessages(require(DescTools))    
suppressMessages(require(data.table))
suppressMessages(require(stringr))
suppressMessages(require(argparser))

#parse arguments:
parser <- arg_parser('Calculate insulation scores of a 3-column long matrix in a_start, b_start, value format at user specified insulation square sizes. Chromosome must be encoded in the filename.', hide.opts = TRUE)

parser <- add_argument(parser, "--cores", type="integer", default=NULL,
                       help="Indicate how many cores should be used for computation. If not set, data.table reads environment variables and uses all ligcal CPUs available")
parser <- add_argument(parser, "--outfile", default=getwd(),
                       help="Indicate path and filename to which output table should be saved")
parser <- add_argument(parser, "--ISmin", short='--ISmin', type="integer", default=100000,  
                       help="Which is the smallest insulation square size in bp for which IS should be calculated")
parser <- add_argument(parser, "--ISmax", short='--ISmax', type="integer", default=1000000, 
                       help="Which is the largest insulation square size in bp for which IS should be calculated")
parser <- add_argument(parser,"--ISsteps", short='--ISsteps', type="integer", default=10, 
                       help="How many different resolutions between ISmin and ISmax should be calculated")
parser <- add_argument(parser, "--input", nargs=Inf,
                       help='Input at least one long matrix in c(\"A_start\", \"B_start\", \"value\") format or use wildcard \"*\" to process multiple chromosomes')



main <- function() {
  args <- parse_args(parser)
  
  #catch error: no matrix supplied
  if(is.na(args$input[1])){
    stop('NO input matrices supplied to function. Please specify at least one matrix in long format for calculation of IS ')
  }
  
  #Print out how many matrices will be converted
  if(length(args$input)==1){
    print('IS will be calculated for 1 matrix')
  }
  else{
    print(paste('IS will be calculated for', length(args$input), 'matrices'))
  }
  
  #set number of threads to be used by data.table:
  setDTthreads(threads=parser$cores)
  
  #initialize a data frame to which mean contacts can be added
  loop_out <- data.table(chrom = character(), viewpoint=integer(), mean_contacts=double(), IS_distance=integer())
  
  #catch error if no writing permission in outfile directory
  if(file.create(args$outfile) == FALSE){
    stop('No writing permission in the --outfile directory')
  }
  
  #loop through all specified files
  for (file in args$input) {
    #Read in the matrix in long format
    matrix_long <- fread(file)
    
    #catch error if matrix has wrong format
    if(dim(matrix_long)[2] > 3){
      stop('matrix has more than 3 columns. input a matrix in A_start, B_start, value format')
    }
    
    #rename or set column names
    setnames(matrix_long, c('A_start', 'B_start', 'count'))
    
    #determine matrix resolution by computing the greatest common divisor of all start posistions
    matrix_res <- DescTools::GCD(unique(matrix_long$A_start))
    matrix_resolution_kb <- paste0(matrix_res/1000, 'kb')

    #add column which indicates how far appart bins A and B are
    matrix_long$contact_distance = matrix_long$B_start - matrix_long$A_start
    #remove contacts of directly neighboring bins and matrix diagonal
    matrix_long_no_diagonal = matrix_long[contact_distance > matrix_res]

    
    #calculate required IS square sizes: defaults to 2 - 20 x matrix resolution in increments of 2 (=10 IS square sizes)
    IS_size_min = args$ISmin
    IS_size_max = args$ISmax
    IS_size_min_to_max = seq(IS_size_min, 
                             IS_size_max,
                             length.out = args$ISsteps)
    
    
    #catch error if matrix resolution does not allow IS calculation:
    if(IS_size_min <= matrix_res){
      stop('--ISmin is smaller or equal to the matrix resolution')
    }
    if(IS_size_max %% matrix_res != 0 |
       IS_size_min %% matrix_res != 0){
      stop('--ISmax or --ISmin is not a multiple of the matrix resolution')
    }
    
    if(sum(IS_size_min_to_max %% matrix_res) != 0){
      stop('Insulation square size specified through --ISsteps does not match the matrix resolution')
    }
    if(IS_size_max <= IS_size_min){
      stop('--ISmax is smaller or equal than --ISmin')
    }
    
    #determine which chromosome is read in from filename
    splitted_filename <- str_split(file, "_")
    chromosome_ID <- splitted_filename[[1]][str_detect(unlist(splitted_filename), '^chr')]
    print(paste('Read in long format matrix of:', chromosome_ID, 'at', matrix_resolution_kb, 'kb resolution', 'from file:', file))
    
    #loop through all IS distances
    for(IS_distance in IS_size_min_to_max){
      #expected_number_of_values = IS_distance/matrix_res * IS_distance/matrix_res
      print(paste('calculating mean contacts within IS square of size:', IS_distance))
      
      # remove everything that is further away than the IS square
      matrix_long_no_diagonal_filtered <- matrix_long_no_diagonal[contact_distance <= IS_distance * 2] 
      
      #how many values are expected in the insulation square of size IS_distance? Make mean count of values in IS square NA if more than half of them are NA
      expected_value_count <- (IS_distance/matrix_res)^2
      
      #loop through all unique start positions
      for(viewpoint in unique(matrix_long_no_diagonal_filtered$A_start)){
        matrix_sub <- matrix_long_no_diagonal_filtered[A_start >= (viewpoint - IS_distance) &
                                                         B_start <= (viewpoint + IS_distance)  ]
        matrix_sub$off_square = dplyr::case_when(matrix_sub$A_start < viewpoint & matrix_sub$B_start <= viewpoint ~TRUE,
                                                 matrix_sub$A_start >= viewpoint & matrix_sub$B_start > viewpoint ~ TRUE,
                                                 TRUE~FALSE)
        matrix_sub <- matrix_sub[off_square == FALSE]
        
        #add row if more than half of the values in the insulation square are there
        if(length(matrix_sub$count) >= expected_value_count/2){
        loop_out <- rbind(loop_out, 
                          list(chrom = chromosome_ID, 
                               viewpoint = as.integer(viewpoint), 
                               mean_contacts =  mean(matrix_sub$count), 
                               IS_distance=as.integer(IS_distance)))
        }
      }
    }
  }
  
  #calculate mean number of contacts per IS square size, per chromosome for IS calculation:
  mean_contacts_perchrom <- loop_out[ , keyby=.(chrom, IS_distance),
                                       .(mean_perchrom=mean(mean_contacts, na.rm=TRUE))]
  IS_table <- merge(loop_out, mean_contacts_perchrom, all.x=TRUE, by=c('chrom', 'IS_distance')) #left join to allow rowwise computation of the IS per chrom and IS_distance
  IS_table[, 'IS' :=  .(log2(mean_contacts) - log2(mean_perchrom))]
  IS_table[, c('mean_contacts', 'mean_perchrom') := NULL]
  
  #reorder columns
  setcolorder(IS_table, c('chrom', 'viewpoint', 'IS_distance', 'IS'))
  #order chromosomes
  setattr(IS_table$chrom, 'levels', str_sort(unique(IS_table$chrom), numeric=TRUE))
  #order by chrom, viewpoint and IS
  IS_table <- IS_table[order(rank(chrom), viewpoint, IS_distance)]
  
  fwrite(x = IS_table, file=args$outfile, sep = '\t', na = 'NA', quote=FALSE)
}
  
main()