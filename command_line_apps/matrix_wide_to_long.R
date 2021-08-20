suppressMessages(require(DescTools))
suppressMessages(require(data.table))
suppressMessages(require(stringr))
suppressMessages(require(argparser))

parser <- arg_parser('Convert a square matrix into a long one', hide.opts = TRUE)

parser <- add_argument(parser, "--cores", type="integer", default=NULL,
                       help="Indicate how many cores should be used for computation. If not set, data.table reads environment variables and uses all ligcal CPUs available")
parser <- add_argument(parser, "--outfile", default=getwd(),
                       help="For processing of multiple files: Indicate path to which output tables should be saved \n 
                       If only a single file is processed, you can sepcify an path + filename")
parser <- add_argument(parser, "--input", nargs=Inf,
                       help='Specify the path to at least one square matrix or use wildcard \"*\" to process multiple chromosomes')


main <- function() {
  args <- parse_args(parser)
  filenames <- args$input
  
  #catch error: no matrix supplied
  if(is.na(filenames[1])){
    stop('NO input matrices supplied to function. Please specify at least one matrix in long format for calculation of IS ')
  }
  
  #Print out how many matrices will be converted
  if(length(filenames)==1){
    print('1 matrix will be converted to long format')
  }
  if(length(filenames) > 1){
    print(paste(length(filenames), 'matrices will be converted to long format'))
  }
  
  #loop through all specified files
  for (filename in filenames) {
    print(paste('read in file:',filename))
    colnames_wide <- colnames(fread(file = filename, nrows = 1))
    mat_wide <- fread(file = filename, colClasses = list(double=2:length(colnames_wide)))
    
    #catch error: input matrix not 'square'
    if(dim(mat_wide)[1]+1 != dim(mat_wide)[2]){
      stop('Input matrix does not meet shape requirements. Please input a matrix that has coordinates (chr:start-end) as first column so that it has 1 more column than rows ')
    }
    
    #reshape matrix and split location identifiers in 3 separate columns
#    colnames_wide <- colnames(mat_wide) now have this part earlier
    mat_long <-melt(mat_wide, id.vars = colnames_wide[1], measure = colnames_wide[-1], na.rm = TRUE)
#    mat_long_noNA <- mat_long[!is.na(mat_long$value)]
    mat_long[, c('A_chrom', 'A_start','A_end') := tstrsplit(get(colnames_wide[1]), ':|-', type.convert=TRUE)]
    mat_long[, c('B_chrom', 'B_start','B_end') := tstrsplit(variable, ':|-', type.convert=TRUE)]

    #catch error: trans chromosomal matrix
    if(length(unique(mat_long$A_chrom)) != 1 |
       length(unique(mat_long$B_chrom)) != 1){
      stop('matrix has trans chromosomal interactions. Reshaping to long format with 3 column output (start, stop, value) would loose information')
    }
    
    #Determine the matrix resolution
    matrix_resolution_A <- DescTools::GCD(unique(mat_long$A_start))
    matrix_resolution_B <- DescTools::GCD(unique(mat_long$B_start))
    
    #catch error: no common devisor of Bins A and B
    if(matrix_resolution_A != matrix_resolution_B){
      stop('Start positions of bin A and bin B have no common devisor')
    }
    matrix_resolution_kb <- paste0(matrix_resolution_A/1000, 'kb')
    
    #Determine which chromosome
    chromosome_ID_A <- unique(mat_long$A_chrom)
    chromosome_ID_B <- unique(mat_long$B_chrom)
    #catch error: more than 1 chromosome in matrix
    if(length(chromosome_ID_A) > 1|
       length(chromosome_ID_B) > 1){
      stop('matrix column or rownames contain more than one chromosome identifer')
    }
    #catch error: not cis interactions
    if(chromosome_ID_A != chromosome_ID_B){
      stop('Chromosome identiiers of rows and columns not matching')
    }
    
    #select A_start, B_start, value
    mat_long_three_column <- mat_long[,c('A_start','B_start', 'value')]
    mat_long_three_column <- mat_long_three_column[!is.na(mat_long_three_column$value)]#drop empty values
    mat_long_three_column <- mat_long_three_column[mat_long_three_column$A_start != mat_long_three_column$B_start] #remove diagonal
    
    #generate filename:
    filename_without_path <- unlist(lapply(str_split(filename, pattern='/'), tail, 1))
    #filename_without_fileending <- unlist(lapply(str_split(filename_without_path, pattern='\\.'), head, 1))
    filename_long_format <- paste(filename_without_path, 'long.tsv.gz', sep="_")
    #filename_chrom_and_res_specificiation <- paste(filename_without_fileending, chromosome_ID_A, matrix_resolution_kb, 'resolution.tsv.gz', sep='_')
    filename_outpath = paste0(args$outfile, filename_long_format)
    
    #write to output directory
    #catch error: no writing permission in output directory
    if(file.create(filename_outpath) != TRUE){
      stop('no writing permission for the output directory')
    }
    fwrite(mat_long_three_column, filename_outpath, sep = '\t')
  }

}

main()
