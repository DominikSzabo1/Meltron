!#/usr/lib64/R

require(DescTools)    
require(data.table)
require(stringr)


main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  output_dir <- args[2]
#  print(output_dir[1])
  filenames <- args[-1:-2]
# print(filenames)
  stopifnot(args[1] =='-o')
  
  #catch error: no matrix supplied
  if(length(args)==0){
    stop('NO input matrices supplied to function. Please specify at least one matrix in long format for calculation of IS ')
  }
  #Print out how many matrices will be converted
  if(length(filenames)==1){
    print('IS will be calculated for 1 matrix')
  }
  else{
    print(paste('IS will be calculated for', length(filenames), 'matrices'))
  }
  
  #initialize a data frame to which mean contacts can be added
  loop_out <- data.table(chrom = character(), viewpoint=integer(), mean_contacts=double(), IS_distance=integer())
  

  #loop through all specified files
  for (file in filenames) {
    #Read in the matrix in long format
    matrix_long <- fread(file)
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
    IS_size_min = matrix_res * 2
    IS_size_max = matrix_res * 20
    IS_size_min_to_max = seq(IS_size_min, 
                             IS_size_max,
                             length.out = 10)
    
    #determine which chromosome is read in from filename
    splitted_filename <- str_split(file, "_")
    chromosome_ID <- splitted_filename[[1]][str_detect(unlist(splitted_filename), 'chr')]
    print(paste('Read in long format matrix of:', chromosome_ID, 'at', matrix_resolution_kb, 'kb resolution', 'from file:', file))
    
    #loop through all IS distances
    for(IS_distance in IS_size_min_to_max){
      #expected_number_of_values = IS_distance/matrix_res * IS_distance/matrix_res
      print(paste('calculating mean contacts within IS square of size:', IS_distance))
      
      # remove everything that is further away than the IS square
      matrix_long_no_diagonal_filtered <- matrix_long_no_diagonal[contact_distance <= IS_distance * 2] 
      
      #loop through all unique start positions
      for(viewpoint in unique(matrix_long_no_diagonal_filtered$A_start)){
        matrix_sub <- matrix_long_no_diagonal_filtered[A_start >= (viewpoint - IS_distance) &
                                                         B_start <= (viewpoint + IS_distance)  ]
        matrix_sub$off_square = dplyr::case_when(matrix_sub$A_start < viewpoint & matrix_sub$B_start <= viewpoint ~TRUE,
                                                 matrix_sub$A_start >= viewpoint & matrix_sub$B_start > viewpoint ~ TRUE,
                                                 TRUE~FALSE)
        matrix_sub <- matrix_sub[off_square == FALSE]
        #      if(length(matrix_sub$count) == expected_number_of_values){
        loop_out <- rbind(loop_out, 
                          list(chrom = chromosome_ID, 
                               viewpoint = as.integer(viewpoint), 
                               mean_contacts =  mean(matrix_sub$count), 
                               IS_distance=as.integer(IS_distance)))
      }
    }
  }
  
  #calculate mean number of contacts per IS square size, per chromosome for IS calculation:
  mean_contacts_perchrom <- loop_out[,
                                     keyby=.(chrom, IS_distance),
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
  
  #generate filename:
#  filename_without_path <- unlist(lapply(str_split(file, pattern='/'), tail, 1))
#  filename_without_fileending <- unlist(lapply(str_split(filename_without_path, pattern='\\.'), head, 1))
#  filename_chrom_and_res_specificiation <- paste(filename_without_fileending, chromosome_ID, matrix_resolution_kb, 'resolution_IS.tsv.gz', sep='_')
  filename_outpath = paste0(output_dir, 'IS.tsv')
  fwrite(IS_table, filename_outpath, sep = '\t', na = 'NA', quote=FALSE)
}
  
main()