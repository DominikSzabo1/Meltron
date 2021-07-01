suppressMessages(require(DescTools))    
suppressMessages(require(data.table))
suppressMessages(require(stringr))
suppressMessages(require(argparser))

#parse arguments:
parser <- arg_parser('Calculates melting scores for supplied genomic intervals by comparing insulation scores of one cell type versus another. Input and reference files need to be in c(\"chrom\", \"viewpoint\", \"IS_distance\" , \"IS\") format', hide.opts = TRUE)

parser <- add_argument(parser, "--cores", type="integer", default=NULL,
                       help="Indicate how many cores should be used for computation. If not set, data.table reads environment variables and uses all ligcal CPUs available")
parser <- add_argument(parser, "--outfile", default=getwd(),
                       help="Indicate path and filename to which output table should be saved")
parser <- add_argument(parser, "--referenceIS", short='-r',
                       help="Indicate path to file that contains insulation scores of reference (usually ESC)")
parser <- add_argument(parser, "--inputIS", short='-i', 
                       help="Indicate path to file that contains insulation scores of input")
parser <- add_argument(parser, "--genomicIntervals", short='-g', 
                       help="Indicate path to file that contains genomic intervals in c(\"chrom\", \"start\", \"end\", \"ID\") format")
parser <- add_argument(parser, "--cutoff_NA", short='-c', default=0.5,
                       help="Indicate percentage (value between 0 and 0.99) of acceptable missing insulation score values. By default, melting scores of genomic intervals that contain more than 50% of missing values will not be computed")

main <- function() {
  args <- parse_args(parser)
  
  #catch error: no matrix supplied
  if(is.na(args$referenceIS[1]) |
     is.na(args$inputIS[1]) ){
    stop('No input or reference insulation score supplied. Please use the -i and -r flags and specify paths to these files')
  }
  
  #catch error: no matrix supplied
  if(args$cutoff_NA >= 1 | 
     args$cutoff_NA < 0 ){
    stop('NA cutoff out of bounds. Needs to be between 0 (no melting scores removed) to 0.99 (at least 99% of insulation scores over genomic intervals are NA)')
  }
  
  #set number of threads to be used by data.table:
  setDTthreads(threads=parser$cores)

  #catch error if no writing permission in outfile directory
  if(file.create(args$outfile) == FALSE){
    stop('No writing permission in the --outfile directory')
  }
  
  #read in files:
  reference_IS <- fread(args$referenceIS)
  input_IS <- fread(args$inputIS)
  genomic_intervals <- fread(args$genomicIntervals)
  
  
  #catch error if matrix has wrong format
  if(dim(input_IS)[2] != 4 | 
     dim(reference_IS)[2] != 4){
    stop('matrix does not have 4 columns. input a matrix in c(\"chrom\", \"viewpoint\", \"IS_distance\" , \"IS\") format')
  }
  
  #set uniform column names
  setnames(reference_IS, c('chrom', 'viewpoint_start', 'IS_distance', 'IS_ref'))
  setnames(input_IS, c('chrom', 'viewpoint_start', 'IS_distance', 'IS_inp'))
  setnames(genomic_intervals, c('chrom', 'start', 'end', 'ID'))
  
  #merge IS tables of reference and input
  IS_table <- merge(input_IS, reference_IS)
  rm(reference_IS, input_IS)
  
  
  #Determine the matrix resolution
  matrix_resolution <- DescTools::GCD(unique(IS_table$viewpoint_start))
  
  #catch error: matrix resolution is lower than 1000 bp 
  if(matrix_resolution < 1000){
    stop('Matrix resolution is below 1000 bp as determined by greatest common divisor of inputIS$viewpoint and referenceIS$viewpoint. Please make sure that the start viewpoint is expressed in multiples of your matrix resolution')
  }
  
  #extend genomic interval coordinates to match matrix_resolution
  genomic_intervals[, 'start_bin' := start %/% matrix_resolution * matrix_resolution ]
  genomic_intervals[, 'end_bin' := (end %/% matrix_resolution * matrix_resolution) + matrix_resolution ]
  
  
  
  #initialize a data frame to which mean contacts can be added
  loop_out <- data.table(chrom = character(), start=integer(), end=integer(), ID = character(), ks_pval_raw=double(), number_compared_values=integer())
  
  for(element in genomic_intervals$ID){
    element_row <- genomic_intervals[ID==element,]
    IS_table_subset <- IS_table[chrom==element_row$chrom &
                                viewpoint_start >= element_row$start_bin &
                                viewpoint_start < element_row$end_bin  ]
    
    ks_pval_row <- ks.test(IS_table_subset$IS_ref,
                      IS_table_subset$IS_inp,
                      alternative = "less")$p.value
    number_compared_values <- length(IS_table_subset$IS_ref)
    
    loop_out <- rbind(loop_out, list(chrom = element_row$chrom,
                                     start = element_row$start_bin,
                                     end = element_row$end_bin,
                                     ID = element_row$ID,
                                     ks_pval_raw = ks_pval_row,
                                     number_compared_values = number_compared_values))
  }
  
  #perform multiple testing correction
  loop_out[, ks_pval_corrected := p.adjust(ks_pval_raw, method = "bonferroni")]
  #calculate melting score from corrected pval
  loop_out[, melting_score:= format(-log10(ks_pval_corrected), digits=1)]
  loop_out[, expected_values := (end - start) / matrix_resolution * length(unique(IS_table$IS_distance))]
  loop_out[, percent_NAs := 1 - (number_compared_values / expected_values) ]
  #recode melting scores to NA if more than given percentage of bins are missing/ would contain NA values
  loop_out[, melting_score  := fcase(percent_NAs > args$cutoff_NA, 'NA', 
                                     percent_NAs <= args$cutoff_NA, melting_score)]
  file_out <- loop_out[, .(chrom, start, end, ID, melting_score, percent_NAs)]
  
  fwrite(file_out, args$outfile, sep='\t', na='NA', quote=FALSE)
}

main()