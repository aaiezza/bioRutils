#!/usr/bin/Rscript

# # # #
# # Just some random R functions that are helpful # #

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

# Logger color constants
logger.colors <- list(
    NORMAL       = list( bg = '',          fg = ''     ),
    STAGE        = list( bg = 'dark grey', fg = ''     ),
    GENE         = list( bg = '',          fg = 'blue' ),
    FILE_PATH    = list( bg = '',          fg = 208    ),
    CONDITION    = list( bg = '',          fg = 2      ),
    IGNORED_COND = list( bg = '',          fg = 1      )
)

logger.getColors <- function() return( names( logger.colors ) )

logger <- function( ..., level = NORMAL, fg, bg )
{
    message <- ...

    if ( missing( fg ) )
        fgl = logger.colors[[level]]$fg
    else fgl = fg

    if ( missing( bg ) )
        bgl = logger.colors[[level]]$bg
    else bgl = bg

    if ( level == 'STAGE' )
        message <- paste( message, ' ' )

    if ( !missing(  ) )

    cat( xtermStyle::style( message, fg = fgl, bg = bgl ) )

    if ( level == 'STAGE' )
        cat( '\n' )
}
