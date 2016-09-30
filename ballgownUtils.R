#!/usr/bin/Rscript

# # # #
# # Just some functions to use with Ballgown # #

suppressMessages( source( '/cvri/Rutils/randomTools.R' ) )
suppressMessages( library( ballgown ) )
# suppressMessages( library( RSkittleBrewer ) )
suppressMessages( library( genefilter ) )
suppressMessages( library( dplyr ) )
suppressMessages( library( devtools ) )

getBallgownData <- function(
    phenoDataFile,
    coVariate = NULL, adjustVars = c( 'population' ),
    sample_pattern = '\\d[ABC]',
    alpha = 5e-2 )
{
    # Retrieve phenotype data file
    ## TODO logging!
    pheno_data <<- read.delim( phenoDataFile )

    if ( is.null( coVariate ) )
    {
        coVariate <- names( pheno_data )[2]
    }

    bg <<- ballgown( dataDir = 'ballgown', samplePattern = sample_pattern, pData = pheno_data )

    # Filter
    bg_filt <<- subset( bg, 'rowVars(texpr(bg)) >1', genomesubset = TRUE )

    results_transcripts <<- stattest( bg_filt,
        feature = 'transcript', covariate = coVariate,
        adjustvars = adjustVars,
        getFC = TRUE, meas = 'FPKM' )

    results_genes <<- stattest( bg_filt,
        feature = 'gene', covariate = coVariate,
        adjustvars = adjustVars,
        getFC = TRUE, meas = 'FPKM' )

    # Associate transcripts with Gene names
    results_transcripts <<- data.frame(
        geneNames = ballgown::geneNames( bg_filt ),
        geneIDs   = ballgown::geneIDs( bg_filt ), results_transcripts )

    results_transcripts <<- arrange( results_transcripts, qval )
    results_genes <<- arrange( results_genes, qval )

    rT_filt <<- subset( results_transcripts, results_transcripts$pval < alpha )
    rG_filt <<- subset( results_genes, results_genes$pval < alpha )

    rT_filt <<- rT_filt[grepl('[^\\.]', rT_filt$geneNames, perl = TRUE),]
}

# Write out genes
writeBGgenes <- function( dir = 'ballgownOutput' )
{
    dir.create( file.path( dir ), showWarnings = FALSE )
    write.Table( results_transcripts, paste0( dir, '/transcript_results.txt' ) )
    write.Table( results_genes, paste0( dir, '/gene_results.txt' ) )
}
