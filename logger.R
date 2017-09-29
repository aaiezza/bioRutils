#!/usr/bin/Rscript

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
    outside,
    prepend          = if(!missing(outside)) outside else logger.format[[level]]$prepend,
    append           = if(!missing(outside)) outside else logger.format[[level]]$append,
    formatted,
    formattedPrepend = if(!missing(formatted)) formatted else logger.format[[level]]$formattedPrepend,
    formattedAppend  = if(!missing(formatted)) formatted else logger.format[[level]]$formattedAppend )
{
    message <- paste( ... )

    if ( !is.null( level ) && level == logger.levels$FILE_PATH )
        message <- normalizePath( message )

    message <- xtermStyle::style( formattedPrepend, message, formattedAppend,
            fg = fg, bg = bg, sep = '' )

    output <- paste0( prepend, message, append )

    if ( print )
        cat( output )
    else return( output )
}

# Print logger levels
print.logger.levels <- function()
{
        null <- lapply( logger.levels, function( l ) {
                logger( l, '\n ---\n', logger('Testing!', level = l, print = FALSE ), '\n ---\n' )
        } )
}

