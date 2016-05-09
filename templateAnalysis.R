#!/usr/bin/Rscript

## Properties
FRM_GENE_FILE <- 'frm_gene_exp.diff.txt'
FAV_COLS <- c( 'gene', 'sample_1', 'sample_2', 'value_1', 'value_2', 'log2(fold_change)', 'p_value' )
GENE_ID_FILTER <- '^((Mir([\\d]|let)+.*)|---|LOC\\d+|(.*Rik.*))$'
HEATMAPS_DIR <- 'heatmaps'
VOLCANO_WIDGET_DIR <- 'interactiveVolcanoPlots'

## Package Dependencies
suppressMessages( library( xtermStyle ) )
suppressMessages( library( xlsx ) )
suppressMessages( library( plotly ) )
suppressMessages( library( plyr ) )
suppressMessages( library( gplots ) )
suppressMessages( library( data.table ) )

## Other Generalized Scripts
source( '/cvri/Rutils/randomTools.R' )
source( '/cvri/Rutils/plotUtils.R' )
source( '/cvri/Rutils/geneSetUtils.R' )
## Prefered options
options( width = 120, warn = -1 )

##
# Import the RNA Seq file into the global variable geneData
#
importGD <- function( file = FRM_GENE_FILE, ... )
{
    # logger( ' Working in', getwd(), '\n## Importing Gene Data',  )
    cat( xtermStyle::style( ' Working in', getwd(), '\n## Importing Gene Data ', bg = 'dark grey' ), '\n' )

    # Import all of the data
    geneData <<- fread( file, header = TRUE, sep = '\t',
        colClasses = c( rep('character', 6), rep('double', 6), 'character' ), data.table = FALSE )

    geneData <<- narrowDownGenes(
        geneData,
        geneFilter = GENE_ID_FILTER,
        geneSetName = paste( 'Gene Set loaded from:',
            xtermStyle::style( normalizePath( file ), fg = 208 ) ), ... )

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

    cat( xtermStyle::style( '## Peek at data ', bg = 'dark grey' ), '\n ',
        xtermStyle::style( normalizePath( filename ), fg = 208 ), '\n', sep='' )
    cat( ' ', xtermStyle::style( sprintf( '%12s', geneData$sample_1[1] ), fg = 2 ), ':',
        paste( head( S2y ), collapse = ', ' ), '\n' )
    cat( ' ', xtermStyle::style( sprintf( '%12s', geneData$sample_2[1] ), fg = 2 ), ':',
        paste( head( S1y ), collapse = ', ' ), '\n' )

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

    cat( xtermStyle::style( ' # Normalize FPKMs ', bg = 'dark grey' ), '\n' )

    fpkms <- matrix( unlist( geneSet[,2:ncol(geneSet)] ), ncol = ncol(geneSet)-1 )
    geneSet$fpkmZ <- t( apply( fpkms, 1, getFPKM ) )
    return( geneSet )
}

##
# Creates a heatmap to a PNG file
#
makeHeatmap <- function( geneSet, file = 'GOI_expression_heatmap.png',
    width = 2200, height = 2200, title = 'GOI Expression',
    heatmapDirectory = HEATMAPS_DIR,
    cexRow = 2.0, cexCol = 5.0, Colv = FALSE, dendrogram = 'row',
    lhei = c(0.05,0.9), labRow = geneSet$gene,
    hclustF = function(d) hclust(d=d, method='ward.D2'), ... )
{
    if ( ncol( conditionsSet ) <= 1 || nrow( geneSet ) <= 1 )
        return()

    cat( xtermStyle::style( '## Creating Heatmap ', bg = 'dark grey' ), ' ',
        xtermStyle::style( nrow(geneSet), fg = 'blue' ), ' Genes being Mapped\n ',
        xtermStyle::style( normalizePath( file ), fg = 208 ), '\n', sep='' )

    dir.create( file.path( heatmapDirectory ) )
    png( file = paste( heatmapDirectory, file, sep = '/'), width = width, height = height )

    par( cex.main = 2.0 )
    # By FPKM
    heatmap.2( matrix( unlist( geneSet$fpkmZ ), ncol = ncol(geneSet$fpkmZ) ),
        main = title, scale = 'row', key = FALSE, lhei = lhei, lwid = c(0.15,0.85),
        Colv = Colv, dendrogram = dendrogram, labRow = labRow,
        labCol = names( geneSet[,c(2:(ncol(geneSet$fpkmZ)+1))] ),
        col = greenred, na.color = 'grey', trace = 'none',
        margin = c( 40, 15 ), cexRow = cexRow, cexCol = cexCol,
        hclust = hclustF, ... )
    suppressMessages( graphics.off() )
}

##
# Print Genes of Interest to a .tsv file
#
writeGOI <- function( geneData, dir = 'goi', fileBody = 'GOI' )
{
    file <- ffn( dir=dir, outputFile = fileBody, ext='.tsv' )
    cat( xtermStyle::style( '## Write GOI to file ', bg = 'dark grey' ), ' ',
        xtermStyle::style( nrow(geneData), fg = 'blue' ), ' Genes being printed\n ',
        xtermStyle::style( normalizePath( file ), fg = 208 ), '\n', sep='' )

    write.Table( unique( geneData$gene ), file = file )
}

massive <- function( geneData, cases = 2 )
{
    cat( xtermStyle::style( '## Aggregate data to Gene Set ', bg = 'dark grey' ), '\n' )

    # Make massive geneSet
    geneSet <- data.frame()
    # This here retreives all of the expression levels
    #  for each gene for each condition regardless of which sample number it is
    #  in the condition comparison
    for ( c in 1:length(conditions) )
    {
        cond <- geneData[geneData$sample_1 == conditions[c],FAV_COLS]
        geneSet <- rbind.fill( geneSet, setNames( cond[,c('gene', 'value_1')], c( 'gene', conditions[c] ) ) )

        cond <- geneData[geneData$sample_2 == conditions[c],FAV_COLS]
        geneSet <- rbind.fill( geneSet, setNames( cond[,c('gene', 'value_2')], c( 'gene', conditions[c] ) ) )
    }

    # Ordering the gense serves no functional purpose, but makes visualizing easier if debugging is needed
    geneSet <- geneSet[order(geneSet$gene),]

    # The row bindings take place very separately, so the set needs to be aggregated
    geneSet <- aggregate( geneSet[, 2:ncol( geneSet )],
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
    cat( sprintf( '%50s: %s\n%50s: %s\n%50s: %s\n',
        paste( 'Genes with expression in at least', cases, 'conditions' ), xtermStyle::style( nrow( geneSet ), fg = 'blue' ),
        'Genes with expression across all conditions',
            xtermStyle::style( nrow( geneSet.uniFPKM ), fg = 'blue' ),
        'Genes missing expression in at least 1 condition',
            xtermStyle::style( nrow( geneSet.naFPKM ), fg = 'blue' ) ), '\n' )

    if ( cases >= length( conditions ) )
    {
        # TODO
        ## geneSet.naFPKM cannot be calculated in this case
    }

    return( list( all = geneSet, naFPKM = geneSet.naFPKM, uniFPKM = geneSet.uniFPKM ) )
}

##
# Extract SDE Genes and place them in sets
#
getSDEGeneSet <- function( gene_set = geneSet$all,
    alpha = 5e-15, altAlpha = 5e-2,
    l2fccU = 1.0, l2fccD = -1.0, ... )
{
    cat( xtermStyle::style( '## Identify significant genes ', bg = 'dark grey' ), '\n' )
    # enriched <<- geneSet[geneSet$gene_id %in% narrowDownGenes(
    #     geneData[geneData$gene_id %in% geneSet$gene_id,],
    #     significant = 'yes', status = 'OK',
    #     alpha = 5e-3, log2FoldChangeCutoff = l2fccU, exprRegulation = 'UP',
    #     geneSetName = 'Enriched Genes' )$gene_id,]
    # depleted <<- geneSet[geneSet$gene_id %in% narrowDownGenes(
    #     geneData[geneData$gene_id %in% geneSet$gene_id,],
    #     significant = 'yes', status = 'OK',
    #     alpha = 5e-3, log2FoldChangeCutoff = l2fccD, exprRegulation = 'DOWN',
    #     geneSetName = 'Depleted Genes' )$gene_id,]

    # sdeGenes <- geneSet.naFPKM[geneSet.naFPKM$gene_id %in% narrowDownGenes(
    #   geneData[geneData$gene_id %in% geneSet.naFPKM$gene_id,],
    #   alpha = 0.0005, l2fccU = l2fccU, l2fccD = l2fccD,
    #   significant = 'yes', status = 'OK' )$gene_id,]

    # geneSet[geneSet$gene_id %in% narrowDownGenes(
    #     geneData[geneData$gene_id %in% geneSet$gene_id,],
    #     alpha = 5e-15, l2fccU = l2fccU, l2fccD = l2fccD,
    #     significant = 'yes', status = 'OK' )$gene_id,]

    # sdeGenes <- geneSet.uniFPKM[geneSet.uniFPKM$gene_id %in% narrowDownGenes(
    #   geneData[geneData$gene_id %in% geneSet.uniFPKM$gene_id,],
    #   alpha = 0.05, l2fccU = l2fccU, l2fccD = l2fccD,
    #   significant = 'yes' )$gene_id,]


    sdeGenes.data <<- narrowDownGenes(
        geneData[geneData$gene %in% gene_set$gene,],
        alpha = alpha, l2fccU = l2fccU, l2fccD = l2fccD,
        significant = 'yes', geneSetName = 'SDE Gene Set' )
    sdeGenes <<- gene_set[gene_set$gene %in% sdeGenes.data$gene,]
    cat( '  ', xtermStyle::style( nrow(sdeGenes), fg = 'blue' ), ' genes present across conditions\n' )

    enriched.data <<- narrowDownGenes(
        geneData[geneData$gene %in% gene_set$gene,],
        significant = 'yes', alpha = altAlpha,
        log2FoldChangeCutoff = l2fccU, exprRegulation = 'UP',
        geneSetName = 'Enriched Genes' )
    enriched.data <<- enriched.data[enriched.data$sample_1 %in% conditions,]
    enriched.data <<- enriched.data[enriched.data$sample_2 %in% conditions,]
    enriched <<- gene_set[gene_set$gene %in% enriched.data$gene,]
    cat( '  ', xtermStyle::style( nrow(enriched), fg = 'blue' ), ' genes present across conditions\n' )

    depleted.data <<- narrowDownGenes(
        geneData[geneData$gene %in% gene_set$gene,],
        significant = 'yes', alpha = altAlpha,
        log2FoldChangeCutoff = l2fccD, exprRegulation = 'DOWN',
        geneSetName = 'Depleted Genes' )
    depleted.data <<- depleted.data[depleted.data$sample_1 %in% conditions,]
    depleted.data <<- depleted.data[depleted.data$sample_2 %in% conditions,]
    depleted <<- gene_set[gene_set$gene %in% depleted.data$gene,]
    cat( '  ', xtermStyle::style( nrow(depleted), fg = 'blue' ), ' genes present across conditions\n' )
}

##
# Parse Conditions
#
parseCondtions <- function( ignoreConditions = c() )
{
    conditionComparisons <<- unique( geneData[,FAV_COLS][,2:3] )
    allConditions <<- unique(as.vector(unlist(conditionComparisons)))

    igCond <- paste( '^(', paste( ignoreConditions, collapse='|' ), ')$', sep='' )
    conditions <<- allConditions[grep( igCond, allConditions, invert = TRUE )]
    conditionsSet <<- t( combn( conditions, 2 ) )

    cat( xtermStyle::style( '## Identifying comparisons \n  ', bg = 'dark grey' ) )
    for ( c in 1:length(conditions) )
        cat( ifelse( c == 1, '', ',' ), xtermStyle::style( conditions[c], fg = 2 ) )

    if ( length(allConditions[grep( igCond, allConditions )]) > 0 )
        cat( ' ~ ignoring', xtermStyle::style( paste(
            allConditions[grep( igCond, allConditions )],
            collapse = ', '), fg = 1 ) )

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

    parseCondtions( ignoreConditions )

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
    else cat( xtermStyle::style( '## Volcano plot comparisons \n', bg = 'dark grey' ) )

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

    # Heatmap

    # makeHeatmap( enriched, file = 'GOI_expression_heatmap_ENRICH.png',
    #   title = 'GOI Significantly Enriched' )
    # makeHeatmap( depleted, file = 'GOI_expression_heatmap_DEPLET.png',
    #   title = 'GOI Significantly Depleted' )
    # makeHeatmap( sdeGenes, colsep = colsep, title = 'SDE Genes' )

    # dataPeek( geneData, len = 20000, cap = 50, file = TRUE )

    # library(GOexpress)
    # BP.5 <- subset_scores( result = )

    ## For now this alternative is handled by HOMER
    # writeGOI( enriched.data, fileBody = 'goi_ENRICH' )
    # writeGOI( depleted.data, fileBody = 'goi_DEPLET' )


    ## In bash
    # /cvri/bin/nohupWrapper.sh eaENRICH/time.txt eaENRICH/output.txt findMotifs.pl goi/goi_ENRICH.tsv $1 eaENRICH/ -depth high -p 2
    # /cvri/bin/nohupWrapper.sh eaDEPLET/time.txt eaDEPLET/output.txt findMotifs.pl goi/goi_DEPLET.tsv $1 eaDEPLET/ -depth high -p 2
    # watchIt() { while :; do c; date; echo; ps -p `pgrep -f findMotifs.pl | head -n 1` -o etime="findMotifs.pl running for:"; echo; tail -n $1 */output.txt; echo; sleep 2; done; }
    # watchIt 20
    # cat( '/cvri/bin/kickOffHOMER' )

    # createExcelFile( geneData )
}

##
# Run HOMER for GSEA and Motif Finding
#   (Totes can't figure out how to run it from R...)
#
runHOMER <- function()
{
    ## In bash where $1 is an organism
    # /cvri/bin/nohupWrapper.sh eaENRICH/time.txt eaENRICH/output.txt findMotifs.pl goi/goi_ENRICH.tsv $1 eaENRICH/ -depth high -p 2
    # /cvri/bin/nohupWrapper.sh eaDEPLET/time.txt eaDEPLET/output.txt findMotifs.pl goi/goi_DEPLET.tsv $1 eaDEPLET/ -depth high -p 2
    # watchIt() { while :; do c; date; echo; ps -p `pgrep -f findMotifs.pl | head -n 1` -o etime="findMotifs.pl running for:"; echo; tail -n $1 */output.txt; echo; sleep 2; done; }
    # watchIt 20
    cat( '/cvri/bin/kickOffHOMER' )
}

##
# Analyze the results of HOMER after it runs
#
analyzeEnrichment <- function(
    analysis = c( 'eaENRICH/biological_process.txt',
                    'eaDEPLET/biological_process.txt',
                    'eaENRICH/kegg.txt',
                    'eaDEPLET/kegg.txt' ),
    n = 5 )
{
    cat( xtermStyle::style( '## Analyze HOMER Enrichment Analysis Output ', bg = 'dark grey' ), '\n' )

    ea <<- list()
    for ( i in 1:length(analysis) )
        ea[[analysis[i]]] <<- read.table( analysis[i], sep = '\t', header = TRUE, fill = TRUE )

    # bp_enriched <<- read.table( 'eaENRICH/biological_process.txt', sep = '\t', header = TRUE, fill = TRUE )
    # bp_depleted <<- read.table( 'eaDEPLET/biological_process.txt', sep = '\t', header = TRUE, fill = TRUE )
    # kg_enriched <<- read.table( 'eaENRICH/kegg.txt', sep = '\t', header = TRUE, fill = TRUE )
    # kg_depleted <<- read.table( 'eaDEPLET/kegg.txt', sep = '\t', header = TRUE, fill = TRUE )

    for ( i in 1:length(ea) )
    {
        readName <- gsub( '.txt', '', gsub( '(/)', '-', names(ea)[i], perl = TRUE ) )

        cond <- ea[[i]]
        cond <- cond[cond$Target.Genes.in.Term > 0,]
        cond <- cond[order(cond$logP),]
        cond <- cond[!duplicated(cond[,2]),]

        cat( ' ===>>> ', xtermStyle::style( readName, fg = 13 ), ' <<<===\n-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --\n' );
        print( cond[1:n,c(2,4,6)] )
        cat( '\n' )

        ea[[i]] <<- cond

        write.Table( cond,
            file = ffn( dir='enrichmentAnalysis', outputFile = readName ),
            col.names = TRUE )
    }
}

##
# Write date to Excel file
#
createExcelFile <- function( geneData, file = 'goi.xlsx', tableName = 'Gene Expression Set' )
{
    # Java imports
    IndexedColors  <- 'org.apache.poi.ss.usermodel.IndexedColors'
    # AreaReference  <- 'org.apache.poi.ss.util.AreaReference'
    # CellReference  <- 'org.apache.poi.ss.util.CellReference'

    wb <- createWorkbook( type = 'xlsx' )

    TABLE_ROWNAMES_STYLE <- CellStyle( wb,
            font = Font( wb, name = 'Verdana', heightInPoints = 12 ),
            alignment = Alignment( wrapText = TRUE, horizontal = "ALIGN_CENTER",
                vertical = "VERTICAL_CENTER" ) )
    TABLE_COLNAMES_STYLE <- CellStyle( wb,
            font = Font( wb, name = 'Calibri', isBold = TRUE ),
            alignment = Alignment( wrapText = TRUE, horizontal = "ALIGN_CENTER",
                vertical = "VERTICAL_CENTER" ) )
    TABLE_CELL_STYLE <- CellStyle( wb,
            font = Font( wb, name = 'Consolas', heightInPoints = 12  ),
            alignment = Alignment( vertical = "VERTICAL_CENTER" ) )
    TABLE_CELL_STYLES <- rep( list( TABLE_CELL_STYLE ), ncol( geneData ) )
    names( TABLE_CELL_STYLES ) <- 1:ncol( geneData )

    sheet <- createSheet( wb, sheetName = 'All GOI' )
    ##  Nifty, but we can use something better
    #.jcall(sheet, 'V','setTabColor', as.integer(12) )
    sheet$setTabColor(J(IndexedColors)$GREEN$index)

    # Lay the data down
    addDataFrame( geneData, sheet = sheet, row.names = FALSE,
        startRow = 1, startColumn = 1,
        colnamesStyle = TABLE_COLNAMES_STYLE,
        rownamesStyle = TABLE_ROWNAMES_STYLE,
        colStyle = TABLE_CELL_STYLES )

    gdCols <- !names( geneData ) %in% c('log2(fold_change)', 'test_stat')
    for( i in 1:length( gdCols ) )
        if( gdCols[i] ) autoSizeColumn( sheet, i )

    addAutoFilter( sheet, paste( 'A1:', LETTERS[ncol(geneData)], nrow(geneData), sep='' ) )

    saveWorkbook( wb, file )
}
