#!/usr/bin/Rscript

# # # #
# # Secreted Proteins database # #

suppressMessages( library( data.table ) )

getSPD <- function( species = c() )
{
    cat( xtermStyle::style( paste( '## Importing Secreted Genes ' ), bg = 'dark grey' ), '\n' )
    spd <- fread( '/cvri/data/secreted_proteins_database/spd.name.nr90.species',
            header = TRUE, sep = '\t', data.table = FALSE )

    if ( length(species) <= 0 )
        species <- unique( spd$species )

    spd <- spd[spd$species %in% species,]
    if ( nrow(spd) <= 0 )
        warning( 'No species called (', paste( species, collapse = ', ' ), ')' )

    cat( ' ', xtermStyle::style( nrow(spd), fg = 'blue' ),
        ' Genes in database of Secreted Genes in (',
        paste( species, collapse = ', ' ), ') ', '\n', sep='' )

    return( spd )
}
