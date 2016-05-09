#!/usr/bin/Rscript

suppressMessages( source( '/cvri/Rutils/randomTools.R' ) )

# # # #
# # Gene set retreival, Filtering and Selection # #

# Read in differentially expressed genes lists before doing this
#  Ex: genes <- read.table(
#        "/home/aaiezza/CVRI/lowenstein/BPCM-123_vs_BPC-123.tsv", sep="\t", header=TRUE )
narrowDownGenes <- function(
    genes, alpha, log2FoldChangeCutoff, l2fccU = 0, l2fccD = 0, significant,
    status, value1 = 0, value2 = 0, geneFilter, exprRegulation = '', geneSetName = 'Gene Set',
    print = TRUE )
{
    # Gotta start with something
    sigGenes <- genes

    ## Filter by p_value (non-SDE genes (non-significantly differentially expressed genes))
    if ( !missing( alpha ) )
        sigGenes <- sigGenes[sigGenes$p_value <= alpha,]

    ## Filter by log2( foldChange ) between conditions
    if ( !missing( log2FoldChangeCutoff ) )
        sigGenes <- switch ( exprRegulation,
            UP = sigGenes[sigGenes$`log2(fold_change)` > log2FoldChangeCutoff,],
            DOWN = sigGenes[sigGenes$`log2(fold_change)` < log2FoldChangeCutoff,],
            sigGenes[abs( sigGenes$`log2(fold_change)` ) > log2FoldChangeCutoff,]
        )

    sigGenes <- sigGenes[
            sigGenes$`log2(fold_change)` > l2fccU |
            sigGenes$`log2(fold_change)` <= l2fccD,]

    ## Filter by significance value
    if ( !missing( significant ) )
        sigGenes <- sigGenes[sigGenes$significant == significant,]

    ## Filter by status
    if ( !missing( status ) )
        sigGenes <- sigGenes[sigGenes$status == status,]

    ### TODO: Should maybe be altered in the future
    ####  to be able to filter both above and below a value.
    ## Filter by value_1
    sigGenes <- sigGenes[sigGenes$value_1 > value1,]

    ### TODO: Should maybe be altered in the future
    ####  to be able to filter both above and below a value.
    ## Filter by value_2
    sigGenes <- sigGenes[sigGenes$value_2 > value2,]

    ## Filter genes with a given regular expression
    if( !missing( geneFilter ) )
        sigGenes <- sigGenes[!grepl( geneFilter, sigGenes$gene, perl = TRUE ),]

    # Show percentage of genes that were significantly differentially expressed
    if ( print )
    {
        logger( prepend = paste( '     ', geneSetName, '\n    [' ),
            logger( sprintf( "%7d", nrow( sigGenes ) ), level = logger.levels$GENE, print = FALSE ), '/',
            logger( sprintf( "%7d", nrow(   genes  ) ), level = logger.levels$GENE, print = FALSE ),
            append = sprintf( '] %.2f%% of genes were preserved\n', ( nrow( sigGenes ) / nrow( genes ) ) * 100 ) )
    }

    return( sigGenes[order(sigGenes$p_value),] )
}

# # # #
# # Different Sortings and Printings # #

printTopFoldChange <- function(
    geneSet, i0 = 0, n = 5, upRegulated = TRUE,
    displayedColumns = c( "gene", "log2(fold_change)" ),
    alsoColumns = c(), displayAllColumns = FALSE
)
{
    foldChange <- ifelse( upRegulated, "UP", "DOWN" )
    description <- paste( "Most tremendously ", foldChange, "-Regulated Genes", sep = "" )

    printGOI( geneSet = geneSet, i0 = i0, n = n, decreasing = upRegulated,
        description = description, sortingCol = "log2(fold_change)",
        displayedColumns = unique( c( displayedColumns, alsoColumns ) ),
        displayAllColumns = displayAllColumns
    )
}

printGOI <- function(
    geneSet, i0 = 0, n = 5, decreasing = FALSE,
    description = "Sorted Gene Set", sortingCol = "gene_id",
    displayedColumns = c( "gene", sortingCol ), displayAllColumns = FALSE
)
{
    for ( cond in names( geneSet ) )
    {
        line <- " -"
        for ( i in 1:18 ){ line <- paste( line, "-" ) }

        format <- paste(
            "\n  ", description, "\n    Condition:  %20s\n",
            line, "\n", sep = ""
        )

        # Alias for the current gene set
        gs <- geneSet[[cond]]

        if( displayAllColumns )
            displayedColumns <- names(gs)

        cat( sprintf( format, cond ) )
        print( tail( head(
            gs[ order(
                gs[sortingCol], decreasing = decreasing ),
                displayedColumns],
            n = n + i0 ), n = n ),
        row.names = FALSE )
    }
}

# Combines a named list of different combined conditions
#  which are themselves named lists but of course of each
#  versused samples' expression data
combineGeneSets <- function( geneSet )
{
    geneData <- data.frame()
    for( cond in names(geneSet) )
    {
        tmp <- data.frame( geneSet[[cond]] )
        tmp$sample_comparison <- cond
        geneData <- rbind( geneData, tmp )
    }

    return( geneData[,unique( c( "sample_comparison", names( geneData ) ) )] )
}
