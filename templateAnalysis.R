#!/usr/bin/Rscript

## Properties
FRM_GENE_FILE <- 'frm_gene_exp.diff.txt'
FAV_COLS <- c( 'gene', 'sample_1', 'sample_2', 'value_1', 'value_2', 'log2(fold_change)', 'q_value' )
GENE_ID_FILTER <- '^((Mir([\\d]|let)+.*)|---|LOC\\d+|(.*Rik.*))$'
HEATMAPS_DIR <- 'heatmaps'
VOLCANO_WIDGET_DIR <- 'interactiveVolcanoPlots'

## Package Dependencies
# suppressMessages( library( xlsx ) )
suppressMessages( library( plotly ) )
suppressMessages( library( plyr ) )
suppressMessages( library( gplots ) )
suppressMessages( library( data.table ) )

## Other Generalized Scripts
source( '/cvri/Rutils/randomTools.R' )
source( '/cvri/Rutils/plotUtils.R' )
source( '/cvri/Rutils/geneSetUtils.R' )
source( '/cvri/Rutils/enrichmentAnalysis.R' )
# source( '/cvri/Rutils/xlsxWriter.R' )
## Prefered options
options( width = 120, warn = -1 )

##
# Import the RNA Seq file into the global variable geneData
#
importGD <- function( file = FRM_GENE_FILE, ... )
{
    logger( 'Working in', getwd(), level = logger.levels$NOTIFY )
    logger( 'Importing Gene Data', level = logger.levels$STAGE  )

    # Import all of the data
    geneData <<- fread( file, header = TRUE, sep = '\t',
        colClasses = c( rep('character', 6), rep('double', 6), 'character' ), data.table = FALSE )

    geneData <<- narrowDownGenes(
        geneData,
        geneFilter = GENE_ID_FILTER,
        geneSetName = logger( prepend = 'Gene Set loaded from:', file, level = logger.levels$FILE_PATH, print = FALSE ), ... )

    print( summary(geneData[,c('value_1', 'value_2','log2(fold_change)','p_value')]) )
}

##
# Peek at our data
#
dataPeek <- function( geneData, len = 100, cap = len, log = '',
    file = FALSE, filename = 'dataPeek.png',
    c1 = geneData$value_1, c2 = geneData$value_2 )
{
    gD <- c( c1, c2 )
    x <- seq( min(gD), ceiling(max(gD)), length = len )
    S1groups <- cut(c1, x)
    S2groups <- cut(c2, x)
    S1y <- tapply(c1, S1groups, length)
    S2y <- tapply(c2, S2groups, length)

    S1Col <- rgb(.2,.2,.8, 0.8)
    S2Col <- rgb(.8,.2,.2, 0.8)

    # I want to bin the expression levels for each sample
    if ( file )
        png( file = filename, width = 1000, height = 880 )
    barplot( rbind( S1y[1:cap], S2y[1:cap] ), xlab ="FPKM", ylab = "frequency",
      main = paste( "FPKM Distribution\n(", cap, ' / ', len, ' bins)', sep = '' ), margins = c(10, 10),
      beside = TRUE, col = c( S1Col, S2Col ), log = log,
      legend = c( as.vector(geneData$sample_1)[1], as.vector(geneData$sample_2)[1]  ) )


    logger( 'Peek at data', level = logger.levels$STAGE )
    logger( filename, append = '\n', level = logger.levels$FILE_PATH )
    logger( prepend = ' ', sprintf( '%12s', geneData$sample_1[1] ), level = logger.levels$CONDITION )
    logger( ':', paste( head( S1y ), collapse = ', ' ), level = logger.levels$NORMAL )
    logger( prepend = ' ', sprintf( '%12s', geneData$sample_2[1] ), level = logger.levels$CONDITION )
    logger( ':', paste( head( S2y ), collapse = ', ' ), level = logger.levels$NORMAL )

    if( file )
        suppressMessages( graphics.off() )
}

##
# Add FPKM Z scores to a geneSet
#
allfpkmZscores <- function( geneSet )
{
    getFPKM <- function( row )
    {
        meanFPKM <- mean( row[complete.cases(row)] )
         sdFPKM  <- sd( row[complete.cases(row)] )

        fpkmRow <- c()
        for( i in 1:length( row ) )
        {
            fpkmRow[i] <- ifelse( is.nan( row[i] ),
                NaN,
                ifelse( is.na( sdFPKM ), 0, ( row[i] - meanFPKM ) / sdFPKM )
            )
        }

        return ( fpkmRow )
    }

    logger( 'Normalize FPKMs', level = logger.levels$STAGE, append = ' ' )
    logger( nrow(geneSet), level = logger.levels$GENE, append = '\n' )

    fpkms <- matrix( unlist( geneSet[,2:ncol(geneSet)] ), ncol = ncol(geneSet)-1 )
    geneSet$fpkmZ <- t( apply( fpkms, 1, getFPKM ) )
    return( geneSet )
}

##
# Creates a heatmap to a PNG file
#
makeHeatmap <- function( geneSet, file = 'GOI_expression_heatmap.png',
    width = 2200, height = 2200, title = 'GOI Expression',
    heatmapDirectory = HEATMAPS_DIR, byFoldChange = FALSE,
    cexRow = 2.0, cexCol = 5.0, Colv = FALSE, dendrogram = 'row',
    lhei = c(0.05,0.9), labRow = geneSet$gene,
    labCol = names( subset( geneSet, select=-c(gene,fpkmZ) ) ),
    hclustF = function(d) hclust(d=d, method='ward.D2'), ... )
{
    if ( ncol( conditionsSet ) <= 1 || nrow( geneSet ) <= 1 )
        return()

    if( !dir.exists( heatmapDirectory ) )
        dir.create( heatmapDirectory, recursive = TRUE )

    logger( 'Creating Heatmap', level = logger.levels$STAGE, append = ' ' )
    logger( nrow(geneSet), level = logger.levels$GENE, append = ' Genes being Mapped\n' )
    logger( file, level = logger.levels$FILE_PATH, append = '\n' )

    dir.create( file.path( heatmapDirectory ) )
    png( file = paste( heatmapDirectory, file, sep = '/'), width = width, height = height )

    par( cex.main = 2.0 )
    # By FPKM
    if ( !byFoldChange && 'fpkmZ' %in% colnames( geneSet )  )
        heatmap.2( as.matrix( geneSet$fpkmZ ),
            main = title, key = FALSE, lhei = lhei, lwid = c(0.15,0.85),
            Colv = Colv, dendrogram = dendrogram, labRow = labRow,
            labCol = labCol,
            col = greenred, na.color = 'grey', trace = 'none',
            margin = c( 40, 15 ), cexRow = cexRow, cexCol = cexCol,
            hclust = hclustF, ... )
    else
        heatmap.2( as.matrix( subset( geneSet, select=-gene ) ),
            main = title, scale = 'row', key = FALSE, lhei = lhei, lwid = c(0.15,0.85),
            Colv = Colv, dendrogram = dendrogram, labRow = labRow,
            col = greenred, na.color = 'grey', trace = 'none',
            margin = c( 40, 15 ), cexRow = cexRow, cexCol = cexCol,
            hclust = hclustF, ... )
    suppressMessages( graphics.off() )
}

##
# Print Genes of Interest to a text-separated-values file
#
writeGOI <- function( geneData, withScore = FALSE, dir = 'goi', fileBody = 'GOI' )
{
    file <- ffn( dir=dir, outputFile = fileBody, ext='.txt' )
    file.create( file )

    logger( 'Write GOI to file', level = logger.levels$STAGE, append = ' ' )
    logger( nrow( geneData ), level = logger.levels$GENE, append = ' Genes being printed\n' )
    logger( file, level = logger.levels$FILE_PATH, append = '\n' )

    if ( 'log2(fold_change)' %in% names( geneData ) )
    {
        if ( withScore )
            write.Table( aggregate( geneData[,c( 'gene', 'log2(fold_change)' )],
                    by = geneData[, 'gene', drop = FALSE], FUN = mean )[,-2],
                file = file, col.names = FALSE, verbose = FALSE )
        else
            write.Table( unique( geneData$gene ), file = file, col.names = FALSE, verbose = FALSE )
    } else
    {
        ## if here, we are assuming that the data given was a matrix of fpkms
        write( theConditions, ncolumns = length( theConditions ), file = file )
        write.Table( geneData[,-length( names( geneData ) )], file = file, append = TRUE, col.names = FALSE, verbose = FALSE )
    }
}

massive <- function( geneData, cases = 2 )
{
    logger( 'Aggregate data to Gene Set', level = logger.levels$STAGE )

    # Make massive geneSet
    geneSet <- data.frame()
    # This here retreives all of the expression levels
    #  for each gene for each condition regardless of which sample number it is
    #  in the condition comparison
    for ( c in 1:length(theConditions) )
    {
        cond <- geneData[geneData$sample_1 == theConditions[c],FAV_COLS]
        geneSet <- rbind.fill( geneSet, setNames( cond[,c('gene', 'value_1')], c( 'gene', theConditions[c] ) ) )

        cond <- geneData[geneData$sample_2 == theConditions[c],FAV_COLS]
        geneSet <- rbind.fill( geneSet, setNames( cond[,c('gene', 'value_2')], c( 'gene', theConditions[c] ) ) )
    }

    # Ordering the genes serves no functional purpose, but makes visualizing easier if debugging is needed
    # geneSet <- geneSet[order(geneSet$gene),]

    # The row bindings take place very separately, so the set needs to be aggregated
    geneSet <- aggregate( geneSet[, -1],
                    by = geneSet[, 'gene', drop = FALSE],
                    mean,  na.rm = TRUE, na.action = NULL )

    # Add the Z scores
    geneSet <- allfpkmZscores( geneSet )
    geneSet <- geneSet[apply(geneSet$fpkmZ, 1, function(x) length( x[complete.cases(x)] ) >= cases ),]

    # Not every condition has expression for every gene.
    #  We may still be interested however in the genes that were only expressed
    #  in less than than all of the conditions.
    #  We can separate them here:
    geneSet.uniFPKM <- geneSet[complete.cases( geneSet ),]
    geneSet.naFPKM  <- geneSet[!( geneSet$gene %in% geneSet.uniFPKM$gene ),]

    logger( sprintf( '%50s: %s\n%50s: %s\n%50s: %s\n',
        paste( 'Genes with expression in at least', cases, 'conditions' ),
            logger( nrow( geneSet ), level = logger.levels$GENE, print = FALSE ),
        'Genes with expression across all conditions',
            logger( nrow( geneSet.uniFPKM ), level = logger.levels$GENE, print = FALSE ),
        'Genes missing expression in at least 1 condition',
            logger( nrow( geneSet.naFPKM ), level = logger.levels$GENE, print = FALSE ) ), append = '\n' )

    if ( cases >= length( theConditions ) )
    {
        # TODO
        ## geneSet.naFPKM cannot be calculated in this case
        logger( 'geneSet.naFPKM cannot be calculated', level = logger.levels$ERROR )
    }

    return( list( all = geneSet, naFPKM = geneSet.naFPKM, uniFPKM = geneSet.uniFPKM ) )
}

##
# Extract SDE Genes and place them in sets
#
getSDEGeneSet <- function( gene_set = geneSet$all,
    alpha = 5e-2, altAlpha = 5e-2,
    l2fccU = 1.0, l2fccD = -1.0, ... )
{
    logger( 'Identify significant genes ', level = logger.levels$STAGE )

    sdeGenes.data <<- narrowDownGenes(
        geneData[geneData$gene %in% gene_set$gene,],
        alpha = alpha, l2fccU = l2fccU, l2fccD = l2fccD,
        significant = 'yes', geneSetName = 'SDE Gene Set' )
    sdeGenes <<- gene_set[gene_set$gene %in% sdeGenes.data$gene,]
    logger( '  ', nrow(sdeGenes), append = ' genes present across conditions\n', level = logger.levels$GENE )

    enriched.data <<- narrowDownGenes(
        geneData[geneData$gene %in% gene_set$gene,],
        significant = 'yes', alpha = altAlpha,
        log2FoldChangeCutoff = l2fccU, exprRegulation = 'UP',
        geneSetName = 'Enriched Genes' )
    enriched.data <<- enriched.data[enriched.data$sample_1 %in% theConditions,]
    enriched.data <<- enriched.data[enriched.data$sample_2 %in% theConditions,]
    enriched <<- gene_set[gene_set$gene %in% enriched.data$gene,]
    logger( '  ', nrow(enriched), append = ' genes present across conditions\n', level = logger.levels$GENE )

    depleted.data <<- narrowDownGenes(
        geneData[geneData$gene %in% gene_set$gene,],
        significant = 'yes', alpha = altAlpha,
        log2FoldChangeCutoff = l2fccD, exprRegulation = 'DOWN',
        geneSetName = 'Depleted Genes' )
    depleted.data <<- depleted.data[depleted.data$sample_1 %in% theConditions,]
    depleted.data <<- depleted.data[depleted.data$sample_2 %in% theConditions,]
    depleted <<- gene_set[gene_set$gene %in% depleted.data$gene,]
    logger( '  ', nrow(depleted), append = ' genes present across conditions\n', level = logger.levels$GENE )
}

##
# Parse Conditions
#
parseConditions <- function( ignoreConditions = c() )
{
    conditionComparisons <<- unique( geneData[,FAV_COLS][,2:3] )
    allConditions <<- unique(as.vector(unlist(conditionComparisons)))

    igCond <- paste( '^(', paste( ignoreConditions, collapse='|' ), ')$', sep='' )
    theConditions <<- allConditions[grep( igCond, allConditions, invert = TRUE )]
    conditionsSet <<- t( combn( theConditions, 2 ) )

    logger( 'Identifying comparisons', level = logger.levels$STAGE )
    for ( c in 1:length(theConditions) )
        cat( ifelse( c == 1, '', ',' ), logger( theConditions[c], level = logger.levels$CONDITION, print = FALSE ) )

    if ( length(allConditions[grep( igCond, allConditions )]) > 0 )
        logger( prepend = ' ~ ignoring ',
            paste( allConditions[grep( igCond, allConditions )], collapse = ', '),
            level = logger.levels$IGNORED_COND )

    cat( '\n' )
}

##
# Handles the importing and identification of condtions
#  in the given gene expression data set
#
prerun <- function( uniformExpression = FALSE,
    ignoreConditions = c(), ... )
{
    importGD()

    parseConditions( ignoreConditions )

    geneSet <<- massive( geneData )

    # Can also be run on only genes that
    #  are expressed across all conditions
    if( uniformExpression )
        getSDEGeneSet( gene_set = geneSet$uniFPKM, ... )
    else getSDEGeneSet( gene_set = geneSet$all, ... )
}

# Volcano Plots
volcanos <- function( toPdf = FALSE, myOwnPdf = FALSE,
    volcanoWidgetDir = VOLCANO_WIDGET_DIR,
    geneLists = list(), condSet = conditionsSet, ... )
{
    plots <- list()

    if ( toPdf )
        startPlot( outputFile = 'volcanoPlot.pdf', ... )
    else logger( 'Volcano plot comparisons', level = logger.levels$STAGE )

    for( i in 1:nrow(condSet) )
    {
        if ( length(geneLists) >= i )
            geneList <- geneLists[[i]]
        else geneList <-  c()

        plot <- produceVolcanoPlot(
            geneData[geneData$sample_1 == condSet[i,1] &
                    geneData$sample_2 == condSet[i,2],],
            geneList = geneList, print = toPdf || myOwnPdf,
            title = paste( condSet[i,1], 'vs', condSet[i,2] ), ... )

        plots[[i]] <- plot
    }

    if ( toPdf || myOwnPdf ) suppressMessages( graphics.off() )

    return( plots )
}

##
# Run HOMER for GSEA and Motif Finding
#   (Totes can't figure out how to run it from R...)
#
runHOMER <- function()
{
    # watchIt findMotifs.pl echo 'tail -20 */output.txt'
    cat( '/cvri/bin/kickOffHOMER' )
}
