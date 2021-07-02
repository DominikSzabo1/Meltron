suppressMessages(require(DescTools))    
suppressMessages(require(data.table))
#suppressMessages(require(stringr))
suppressMessages(require(argparser))

#parse arguments:
parser <- arg_parser('Calculates melting scores for supplied genomic intervals by comparing insulation scores of one cell type versus another. Input and reference files need to be in c(\"chrom\", \"viewpoint\", \"IS_distance\" , \"IS\") format', hide.opts = TRUE)

parser <- add_argument(parser, "--cores", type="integer", default=NULL,
                       help="Indicate how many cores should be used for computation. If not set, data.table reads environment variables and uses all ligcal CPUs available")
parser <- add_argument(parser, "--outfile",  default=paste0(getwd(), '/melting_scores.tsv'),
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
    stop('No input or reference insulation score supplied. Please use the -i and -r flags and specify paths to these files or run \"Rscript Meltron.R --help\" to get help')
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
  IS_table <- merge(input_IS, reference_IS, all=TRUE)
  rm(reference_IS, input_IS)
  
  
  #Determine the matrix resolution and IS length scales
  matrix_resolution <- DescTools::GCD(unique(IS_table$viewpoint_start))
  IS_sizes <- unique(IS_table$IS_distance)
  n_IS_sizes <- length(IS_sizes)
  
  #catch error: matrix resolution is lower than 1000 bp 
  if(matrix_resolution < 1000){
    stop('Matrix resolution is below 1000 bp as determined by greatest common divisor of inputIS$viewpoint and referenceIS$viewpoint. Please make sure that the start viewpoint is expressed in multiples of your matrix resolution')
  }
  
  #extend genomic interval coordinates to match matrix_resolution
  genomic_intervals[, 'start_bin' := start %/% matrix_resolution * matrix_resolution ]
  genomic_intervals[, 'end_bin' := (end %/% matrix_resolution * matrix_resolution) + matrix_resolution ]
  
  
  
  #initialize a data frame to which mean contacts can be added
  loop_out <- data.table(chrom = character(), start=integer(), end=integer(), ID = character(), ks_pval_raw=double(), perc_na_or_missing_ref=integer(), perc_na_or_missing_inp=integer())
  
  for(element in genomic_intervals$ID){
    element_row <- genomic_intervals[ID==element,]
    IS_table_subset <- IS_table[chrom==element_row$chrom &
                                viewpoint_start >= element_row$start_bin &
                                viewpoint_start < element_row$end_bin  ]
    n_expected_bins <- (element_row$end_bin - element_row$start_bin) / matrix_resolution
    n_expected_values <- n_expected_bins * n_IS_sizes
    #when NA in both ref and inp, these entries would be missing from the table
    missing <- n_expected_values - dim(IS_table_subset)[1] 
    na_ref <- sum(is.na(IS_table_subset$IS_ref))
    na_inp <- sum(is.na(IS_table_subset$IS_inp))
    perc_na_or_missing_ref <- (missing + na_ref) / n_expected_values
    perc_na_or_missing_inp <- (missing + na_inp) / n_expected_values
    
    #only run KS test if more than x percent NAs in the genomic interval. rows are
    if( args$cutoff_NA  >=  perc_na_or_missing_ref &
        args$cutoff_NA  >=  perc_na_or_missing_inp){
      
      ks_pval_row <- ks.test(IS_table_subset$IS_ref,
                             IS_table_subset$IS_inp,
                             alternative = "less")$p.value
#      number_compared_values_ref <- length(IS_table_subset$IS_ref)
#     number_compared_values_inp <- length(IS_table_subset$IS_inp)

      loop_out <- rbind(loop_out, list(chrom = element_row$chrom,
                                     start = element_row$start_bin,
                                     end = element_row$end_bin,
                                     ID = element_row$ID,
                                     ks_pval_raw = ks_pval_row,
                                     perc_na_or_missing_ref = sprintf(perc_na_or_missing_ref, fmt='%#.2f'),
                                     perc_na_or_missing_inp = sprintf(perc_na_or_missing_inp, fmt='%#.2f')))
    }
    #if too many values missing => make pvalue NA
    else{
      loop_out <- rbind(loop_out, list(chrom = element_row$chrom,
                                       start = element_row$start_bin,
                                       end = element_row$end_bin,
                                       ID = element_row$ID,
                                       ks_pval_raw = NA,
                                       perc_na_or_missing_ref = sprintf(perc_na_or_missing_ref, fmt='%#.2f'),
                                       perc_na_or_missing_inp = sprintf(perc_na_or_missing_inp, fmt='%#.2f')))
    }
  }
  
  #perform multiple testing correction
  loop_out[, ks_pval_corrected := p.adjust(ks_pval_raw, method = "bonferroni")]
  #calculate melting score from corrected pval
  loop_out[, melting_score:= sprintf(-log10(ks_pval_corrected), fmt='%#.2f')]
  file_out <- loop_out[, .(chrom, start, end, ID, melting_score, perc_na_or_missing_ref, perc_na_or_missing_inp)]
  
  fwrite(file_out, args$outfile, sep='\t', na='NA', quote=FALSE)
}

main()