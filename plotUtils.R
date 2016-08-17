#!/usr/bin/Rscript

suppressMessages( library( calibrate ) )
suppressMessages( library( ggplot2 ) )
suppressMessages( library( ggrepel ) )
suppressMessages( library( ggdendro ) )
suppressMessages( library( dendextend ) )
suppressMessages( require( htmlwidgets ) )
suppressMessages( require( gridExtra ) )

Sys.setenv("plotly_username" = "aaiezza")
Sys.setenv("plotly_api_key" = "7igknczkw2")

suppressMessages( source( '/cvri/Rutils/randomTools.R' ) )

# # # #
# # Shortcut for not needing to think to hard on a plot file name

startPlot <- function(
    outputFile = format( Sys.time(), "%Y-%m-%d_%H%M%S_outputPlot.pdf" ),
    dir = '.', width = 15, height = 15, pointsize = 16,
    landscape = FALSE, onefile = TRUE, verbose = TRUE, ... )
{
    if ( verbose )
    {
        logger( 'Creating volcano plot file', level = logger.levels$STAGE )
        logger( outputFile, level = logger.levels$FILE_PATH, append = '\n' )
    }

    if ( !missing( dir ) )
    {
        if ( !dir.exists( dir ) )
            dir.create( dir )
        outputFile <- paste0( dir, '/', outputFile )
    }
    pdf( file = outputFile, width = ifelse(landscape, height, width),
        height = ifelse(landscape, height, width), pointsize = pointsize,
        onefile = onefile )
}

# For multiple plots on the same page:
#  par( mfrow = c(2,3) )

produceVolcanoPlot <- function(
    data, title = "Plot", alpha = 5e-2, log2FoldChangeCutoff = 1.0,
    geneList = c(), namingQValueCutoff = alpha, namingLog2FoldChangeCutoff = log2FoldChangeCutoff,
    namingLog2FoldChangeCutoffDown = -namingLog2FoldChangeCutoff,
    namingLog2FoldChangeCutoffUp = namingLog2FoldChangeCutoff,
    alphaFill = 0.9, nudge = 0.2, showLegend = TRUE, scaleOverMedian, plotStuff,
    print = FALSE,
    toWidget = FALSE, volcanoWidgetDir = 'interactiveVolcanoPlots', ...
)
{
    conds <- strsplit( title, ' vs ' )[[1]]
    if ( length( conds ) == 2 )
        logger( prepend = '  ~',
            sprintf( '%20s vs %-20s', conds[1], conds[2] ),
            append = '', fg = 'green' )
    else logger( prepend = '  ~', sprintf( '%30s', title ), append = '', fg = 'green' )

    data$class <- with( data,
        ifelse( abs( `log2(fold_change)` ) < log2FoldChangeCutoff & p_value > alpha,
            'Not Significant',
        ifelse( abs( `log2(fold_change)` ) > log2FoldChangeCutoff & p_value > alpha,
            paste('Not Significant\n  Considerable Log2FoldChange >', log2FoldChangeCutoff),
        ifelse( abs( `log2(fold_change)` ) < log2FoldChangeCutoff & p_value <= alpha,
            paste('Significant <', alpha, '\n  Inadequate Log2FoldChange'),
        ifelse( abs( `log2(fold_change)` ) > log2FoldChangeCutoff & p_value <= alpha &
            significant == 'no',
            'Significant p-value',
        'Significant q-value') ) ) )
    )

    GENE_NAMES = paste( '^(', paste( geneList, collapse = '|' ), ')$', sep = '' )

    labels <- subset( data, grepl( GENE_NAMES, gene ) |
        ( -log10( q_value ) > -log10( namingQValueCutoff ) &
        ( `log2(fold_change)` > namingLog2FoldChangeCutoffUp |
          `log2(fold_change)` < namingLog2FoldChangeCutoffDown )
        & class == 'Significant q-value' ) )

    if ( !missing(scaleOverMedian) )
        scale <- round( abs( median(data$`log2(fold_change)`) ) + scaleOverMedian, 1 )
        if ( is.na( scale ) ) scale = 1

    plot <- ggplot( data, aes( x = `log2(fold_change)`, y = -log10( q_value ), gene = gene ) ) +
      ####
        geom_point( aes( fill = class ), show.legend = showLegend,
            color = 'black', shape = 21, size = 5, stroke = 1.3 ) +
      ####
        scale_fill_manual( name = NULL,
            breaks = c(
                   'Significant q-value',
                   'Significant p-value',
            paste( 'Significant <', alpha, '\n  Inadequate Log2FoldChange'),
            paste( 'Not Significant\n  Considerable Log2FoldChange >', log2FoldChangeCutoff),
                   'Not Significant' ),
            values = c(
            alpha( 'darkgrey', alphaFill ),
            alpha( 'orange'  , alphaFill ),
            alpha( 'red'     , alphaFill ),
            alpha( 'yellow'  , alphaFill ),
            alpha( 'green'   , alphaFill ) ) ) +
      ####
        theme_bw( base_size = 26 ) +
      ####
        ggtitle( title ) +
      ####
        theme(
            legend.key.height = unit( 2, 'lines' ),
            legend.position = 'bottom',
            legend.direction = 'horizontal',
            plot.margin = unit(c(2,2,2,2), 'lines'),
            plot.title = element_text( size=24, vjust=0.5, margin=margin(10,0,12,0) ),
            axis.title = element_text( vjust=0.5, margin=margin(50,50,50,50) ) ) +
      ####
        guides( fill = guide_legend( ncol = 2, byrow = TRUE ) )


    if ( !missing( plotStuff ) ) plot <- plot + plotStuff
    if ( !missing( scaleOverMedian ) ) plot <- plot + xlim( -scale, scale )
    if ( nrow(labels) > 0 )
    {
        plot <- plot + geom_label_repel(
                data = labels,
                aes( label = gene, fill = class ),
                fontface = 'italic', color = 'black',
                show.legend = FALSE,
                size = 6, force = 2,
                label.size = 1.15,
                segment.size = 0.8,
                arrow = arrow(25, unit(0.01,'npc')),
                box.padding = unit( 0.02, 'npc' ),
                point.padding = unit( 0.6, 'lines' ),
                nudge_x = nudge
            )
    }

    if ( print ) print( plot )

    if ( toWidget )
    {
        cat( ' ~ widget' )
        dir.create( file.path( volcanoWidgetDir ) )
        fileName <- paste( gsub(' ', '_', title), '.html', sep='' )

        widget <- as.widget( ggplotly( plot, tooltip = c( 'gene', 'x', 'y' ) ) )

        htmlwidgets::saveWidget( widget, fileName )
        file.rename( fileName, file.path( volcanoWidgetDir, fileName ) )
    }

    cat( '\n' )

    return( plot )
}

produceDendrogram <- function(
    sampleCluster,
    fileName = ffn( prepend = 'dendrogram', ext = '.png' ),
    plot_labs = NULL, k_clusters = 1, ymin = -0.15 )
{
    png( fileName, width = 1080, height = 850 )

    dend <- sampleCluster %>% as.dendrogram %>%
        set( 'branches_k_color', k = k_clusters ) %>%
        set( 'labels_cex', 1.2 )

    plot <- as.ggdend( dend )

    plot <- ggplot( plot ) +
        theme_bw( base_size = 20 ) +
        theme(
            plot.margin = unit(c(2,2,2,2), 'lines'),
            plot.title = element_text( size=24, vjust=0.5, margin=margin(10,0,25,0) ),
            axis.title.y = element_text( vjust=0.5, margin=margin(0,15,0,5) ),
            axis.title.x = element_text( vjust=0.5, margin=margin(15,0,5,0) ),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank() ) +
        ylim( ymin, max( get_branches_heights( dend ) ) ) +
        plot_labs

    print( plot )
    suppressMessages( graphics.off() )
}
