#!/usr/bin/Rscript

# # # #
# # Just some random R functions that are helpful # #

# Override defaults of functions
gsub <- function( pattern, replacement, x, ignore.case = FALSE, perl = TRUE, fixed = FALSE, useBytes = FALSE )
{
    return( gsub(pattern, replacement, x, ignore.case = ignore.case, perl = perl, fixed = fixed, useBytes = useBytes) );
}

# Generate random nucleotide sequence
randNucleotides <- function(
    n, nucleotides = c('A', 'C', 'G', 'T'),
    weights = c(0.25, 0.25, 0.25, 0.25), GCrich = FALSE,
    counts = FALSE, print = FALSE, sequenceId = 'Sequence' )
{
    sequence <- noquote(
            sample( x = nucleotides, size = n, replace = TRUE, prob = weights ) )

    # Report counts
    if ( print ) print( paste( sequence, collapse=", " ) )

    seqFreq <- table( sequence )
    if ( counts ) print( seqFreq )

    if ( GCrich ) cat( sprintf( " GC content : %.2f%%\n", sum( seqFreq[c( 'C', 'G' )] ) * 100 / n ) )

    return( list( id = sequenceId, seq = sequence ) )
}

# Sequence summary
seqSummary <- function( sequence )
{
    return( setNames( count( sequence$seq ), c( 'nucleotide', 'frequency' ) ) )
}

# Change how wide the R console's prints are
wideScreen <- function( howWide = Sys.getenv( "COLUMNS" ) )
{
  options( width = as.integer( howWide ) )
}

# Favorite file name!
ffn <- function(
    dir, prepend, append, ext = '.txt', createParentPath = TRUE,
    outputFile = format( Sys.time(), "%Y-%0m-%0d_%0H%0M%0S" ) )
{
    if ( !missing( prepend ) )
        outputFile <- paste( prepend, outputFile, sep = '_' )

    if ( !missing( dir ) )
    {
        if ( createParentPath && !dir.exists( dir ) )
            dir.create( dir, recursive = TRUE )
        outputFile <- paste( dir, outputFile, sep = '/' )
    }

    if ( !missing( append ) )
        outputFile <- paste( outputFile, append, sep = '_' )

    return( paste( outputFile, ext, sep = '' ) )
}

# I tend to prefer these defaults on my calls to write.table
write.Table <- function(
    x, file = "", append = FALSE, quote = FALSE, sep = "\t",
    eol = "\n", na = "NA", dec = ".", row.names = FALSE,
    col.names = FALSE, qmethod = c("escape", "double"),
    fileEncoding = "" )
{
    return( write.table(
        x = x, file = file, append = append, quote = quote, sep = sep,
        eol = eol, na = na, dec = dec, row.names = row.names,
        col.names = col.names, qmethod = qmethod,
        fileEncoding = fileEncoding ) )
}

## Logger
# Uses xtermstyle
suppressMessages( require( xtermStyle ) )

# Logger color constants
logger.format <- list(
    NORMAL       = list( bg = '',          fg = '',      formattedPrepend = '',     formattedAppend = '',    prepend = '', append = '\n' ),
    STAGE        = list( bg = 'dark grey', fg = '',      formattedPrepend = '## ',  formattedAppend = ' ',   prepend = '', append = '\n' ),
    NOTIFY       = list( bg = 'dark blue', fg = 'white', formattedPrepend = '    ', formattedAppend = ' ',   prepend = '', append = '\n' ),
    GENE         = list( bg = '',          fg = 'blue',  formattedPrepend = '',     formattedAppend = '',    prepend = '', append = ''   ),
    FILE_PATH    = list( bg = '',          fg = 208,     formattedPrepend = ' ',    formattedAppend = '',    prepend = '', append = ''   ),
    CONDITION    = list( bg = '',          fg = 2,       formattedPrepend = '',     formattedAppend = '',    prepend = '', append = ''   ),
    IGNORED_COND = list( bg = '',          fg = 1,       formattedPrepend = '',     formattedAppend = '',    prepend = '', append = ''   ),
    ERROR        = list( bg =  1,          fg = 'white', formattedPrepend = ' ==',  formattedAppend = '== ', prepend = '', append = '\n' )
)

logger.levels <- setNames( as.list( names( logger.format ) ), names( logger.format ) )

logger <- function( ..., level = logger.levels$NORMAL, print = TRUE,
    fg               = logger.format[[level]]$fg,
    bg               = logger.format[[level]]$bg,
    prepend          = logger.format[[level]]$prepend,
    append           = logger.format[[level]]$append,
    formattedPrepend = logger.format[[level]]$formattedPrepend,
    formattedAppend  = logger.format[[level]]$formattedAppend )
{
    message <- paste( ... )

    if ( level == logger.levels$FILE_PATH )
        message <- normalizePath( message )

    message <- xtermStyle::style( formattedPrepend, message, formattedAppend,
            fg = fg, bg = bg, sep = '' )

    output <- paste( prepend, message, append, sep = '' )

    if ( print )
        cat( output )
    else return( output )
}
