suppressMessages(require(DescTools))    
suppressMessages(require(data.table))
suppressMessages(require(stringr))
suppressMessages(require(argparser))

#parse arguments:
parser <- arg_parser('Calculates melting scores for supplied genomic intervals by comparing insulation scores of one cell type versus another.', hide.opts = TRUE)

parser <- add_argument(parser, "--cores", type="integer", default=NULL,
                       help="Indicate how many cores should be used for computation. If not set, data.table reads environment variables and uses all ligcal CPUs available")
parser <- add_argument(parser, "--outfile", default=getwd(),
                       help="Indicate path and filename to which output table should be saved")
parser <- add_argument(parser, "--referenceIS", short='-r',
                       help="Indicate path to file that contains insulation scores of reference (usually ESC)")
parser <- add_argument(parser, "--inputIS", short='-i', 
                       help="Indicate path to file that contains insulation scores of input")



