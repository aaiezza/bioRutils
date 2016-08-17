#!/usr/bin/Rscript

suppressMessages( source( '/cvri/Rutils/randomTools.R' ) )

##
# Analyze the results of HOMER after it runs
#
analyzeHomerEnrichment <- function(
    dir = 'goi', toFiles = FALSE,
    n = 5 )
{
    logger( 'Analyze HOMER Enrichment Analysis Output ', level = logger.levels$STAGE )

    ANALYSIS <-
        list( gobp = 'biological_process.txt',
              gocc = 'cellular_component.txt',
              gomf = 'molecular_function.txt',
              kegg = 'kegg.txt',
              reac = 'reactome.txt',
              wiki = 'wikipathways.txt',
              pfam = 'pfam.txt' )
    HomerEnrichmentAnalysis <- list()

    orderUp <- function( table )
        return( table[order( table$logP ), !( names(table) %in% c( 'Enrichment', 'Genes in Term', 'Target Genes in Term' ) )] )

    if ( toFiles )
        logger( 'HOMER Enrichment Analysis - \n', logger( dirname(dir), level = logger.levels$FILE_PATH, print = FALSE, formattedPrepend = '' ),
                level = logger.levels$STAGE, append  = ':\n ' )

    for( ea in names( ANALYSIS ) )
    {
        eaFile <- normalizePath( paste( dir, ANALYSIS[[ea]], sep = '/') )
        if( !file.exists( eaFile ) ) next

        HomerEnrichmentAnalysis[[ea]] <-
            read.delim( eaFile )

        HomerEnrichmentAnalysis[[ea]] <- orderUp( HomerEnrichmentAnalysis[[ea]] )

        if ( toFiles )
        {
            logger( basename( ffn( dir = dir, prepend = 'tophits', outputFile = ea ) ), level = logger.levels$FILE_PATH, append = ' ' )
            write.Table( HomerEnrichmentAnalysis[[ea]], ffn( dir = dir, prepend = 'tophits', outputFile = ea ), verbose = FALSE )
        }
    }
    cat( '\n' )

    return( HomerEnrichmentAnalysis )
}

##
# Analyze Enrichment Analysis from DAVID
#
analyzeDAVIDEnrichment <- function(
    gene_ids, listName,
    dir = 'goi', toFiles = FALSE,
    email = 'Alessandro_Aiezza@URMC.Rochester.edu',
    url = 'https://david.ncifcrf.gov/webservice/services/DAVIDWebService' )
{
    suppressMessages( require( RDAVIDWebService ) )
    logger( 'DAVID Enrichment Analysis', level = logger.levels$STAGE )

    david <<- DAVIDWebService$new( email = email, url = url )
    good_gene_ids <- gsub( pattern = '\\.\\d+', replacement = '', gene_ids, perl = TRUE )

    data( good_gene_ids )
    # logger( david$addList( inputIds = genes, idType = 'OFFICIAL_GENE_SYMBOL',
    #         listName = listName, listType = 'Gene' ),
    #     level = logger.levels$NORMAL )
    logger( sprintf( '  %.2f%% of given genes were found', david$addList( inputIds = good_gene_ids, idType = 'ENSEMBL_GENE_ID',
            listName = listName, listType = 'Gene' )$inDavid * 100 ),
        level = logger.levels$NORMAL )
    david$setAnnotationCategories( c( 'GOTERM_BP_ALL' ) )#, 'GOTERM_MF_ALL', 'GOTERM_CC_ALL' ) )

    clusterReport <- david$getClusterReport( type = 'Term' )
    goHits <- david$getFunctionalAnnotationChart()

    fav_cols <- -c(1,6:13)
    goHits <- goHits[order(goHits$PValue),fav_cols]

    if ( toFiles )
    {
        logger( dirname(dir), level = logger.levels$FILE_PATH, append = ' :\n ', bg = logger.format$STAGE$bg )

        # Print the goHits
        logger( basename( ffn( dir = dir, outputFile = 'go_terms' ) ), level = logger.levels$FILE_PATH, append = ' ' )
        write.Table( goHits, file = ffn( dir = dir, outputFile = 'go_terms' ), verbose = FALSE )

        # Print the DAGs
        startPlot( dir = dir, outputFile = 'gobp_graphs.pdf', verbose = FALSE )
        logger( 'gobp_graphs.pdf', level = logger.levels$FILE_PATH, append = ' ' )
        for ( clust in  members( clusterReport ) )
        {
            davidGODag <- DAVIDGODag( clust, pvalueCutoff = 5e-2, 'BP' )
            plotGOTermGraph( g = goDag( davidGODag ), r = davidGODag, max.nchar = 45, node.shape = 'ellipse' )
        }
        suppressMessages( graphics.off() )
        cat( '\n' )
    }

    return( goHits )
}

##
# Analyze Enrichment Analysis from GeneTrail 2
#
analyzeGeneTrail2Enrichment <- function(
    dir = 'enrichment', toFiles = FALSE, m = 10 )
{
    # TODO, finish this class and use IT instead of
    #  relying on the user to run the analysis manually on the website
    # gtClient <<- GeneTrail2Client$new()

    ANALYSIS <-
        list( gobp = 'GO_-_Biological_Process.txt',
              gocc = 'GO_-_Cellular_Component.txt',
              gomf = 'GO_-_Molecular_Function.txt',
              kegg = 'KEGG_-_Pathways.txt',
              reac = 'Reactome_-_Pathways.txt',
              wiki = 'WikiPathways.txt' )
    GTEnrichmentAnalysis <- list()

    orderUp <- function( table )
        return( table[order( table$P_value ), !( names(table) %in% c( 'Reference' ) )] )

    if ( toFiles )
        logger( 'GeneTrail2 Enrichment Analysis - \n', logger( dir, level = logger.levels$FILE_PATH, print = FALSE, formattedPrepend = '' ),
                level = logger.levels$STAGE, append  = ':\n' )

    # After downloading and unzipping into the goi folders from GeneTrail2
    for( ea in names( ANALYSIS ) )
    {
        eaFile <- normalizePath( paste( dir, ANALYSIS[[ea]], sep = '/') )
        if( !file.exists( eaFile ) ) next

        GTEnrichmentAnalysis[[ea]] <- read.delim( eaFile )

        names( GTEnrichmentAnalysis[[ea]] ) <- str_replace( names( GTEnrichmentAnalysis[[ea]] ), '\\.', '_' )

        GTEnrichmentAnalysis[[ea]] <- orderUp( GTEnrichmentAnalysis[[ea]] )

        if ( toFiles )
        {
            logger( ' ', basename( ffn( dir = dir, outputFile = ea ) ), level = logger.levels$FILE_PATH, append = ' ' )
            write.Table( GTEnrichmentAnalysis[[ea]], ffn( dir = dir, outputFile = ea ), verbose = FALSE )
        }
    }
    cat( '\n' )

    return( GTEnrichmentAnalysis )

}

suppressMessages( require( RCurl ) )
suppressMessages( require( jsonlite ) )
suppressMessages( require( R6 ) )
##
# Gene Trail 2 Client class
#
GeneTrail2Client <- R6Class(
    'GeneTrail2Client',
    public = list(
        initialize = function(
            username = 'aaiezza',
            password = 'gtpassword',
            sessionId )
        {
            if ( !missing( username ) && !missing( password ) )
                self$login( username, password )
            if ( !missing( sessionId ) )
                private$sessionId <- sessionId
            else private$sessionId <-self$getSessionId()
            cat( private$sessionId, '\n' )
        },
        getAPIdoc = function()
        {
            return( self$doGet( '/api-docs' ) )
        },
        login = function( username, password )
        {
            return( self$doPost( '/user/login?username=%s&password=%s', username, password ) )
        },
        getSessionId = function()
        {
            if ( is.na( private$sessionId ) )
                private$sessionId <- self$doGet( '/session' )[['session']]
            return( private$sessionId )
        },
        uploadGenes = function( genesFile, displayName )
        {
            cat( paste0( 'https://genetrail2.bioinf.uni-sb.de/results.html?session=', private$sessionId, '&show_all_results=true\n' ) )
            return( fromJSON( RCurl::postForm(
                sprintf( paste0( private$GT_URL, '/upload?session=%s' ),
                    private$sessionId ),
                file = fileUpload( genesFile ) ) ) )#, displayName = displayName ) ) )
        },
        getJobAlgorithms = function()
        {
            return( self$doGet( '/job/algorithms' ) )
        },
        getJobAlgorithmParameters = function( algorithm )
        {
            return( self$doGet( '/job/parameters/%s', algorithm ) )
        },
        startJob = function( contact )
        {
            if ( missing( contact ) )
                return( self$doGet( '/job/start?session=%s', private$sessionId ) )
            else
                return( self$doGet( '/job/start?session=%s&contact=%s', private$sessionId, contact ) )
        },
        setupJob = function( algorithm, ... )
        {
            return( self$doPostForm(
                sprintf( paste0( private$GT_URL, '/job/setup/%s?session=%s' ),
                    algorithm, private$sessionId ), ... ) )
        },
        stopJob = function()
        {
            return( self$doGet( '/job/stop?session=%s', private$sessionId ) )
        },
        runGSEA = function( genes, jobName )
        {
            job <- self$uploadGenes( genes, jobName )
            return( job )
        },
        getGT_URL = function()
        {
            return( private$GT_URL )
        },
        doGet = function( endpoint, ... )
        {
            return( fromJSON( RCurl::httpGET(
                sprintf( paste0( private$GT_URL, endpoint ), ... ),
                httpheader = list(
                    Accept         = 'application/json',
                    `Content-type` = 'application/x-www-form-urlencoded' )
                ) ) )
        },
        doPost = function( endpoint, ... )
        {
            return( fromJSON( RCurl::httpPOST(
                sprintf( paste0( private$GT_URL, endpoint ), ... ),
                httpheader = list(
                    Accept         = 'application/json',
                    `Content-type` = 'application/x-www-form-urlencoded' )
                ) ) )
        },
        doPostForm = function( url, ... )
        {
            return( fromJSON( RCurl::postForm( url, style = 'POST', ... ) ) )
        },
        doPut = function( endpoint, content, ... )
        {
            return( fromJSON( RCurl::httpPUT(
                sprintf( paste0( private$GT_URL, endpoint ), ... ),
                content,
                httpheader = list(
                    Accept         = 'application/json',
                    `Content-type` = 'application/x-www-form-urlencoded' )
                ) ) )
        },
        api = list(
            API_DOC = '/api-docs',
            USER    = list(
                login = '/user/login?username=%s&password=%s'
            ),
            SESSION = list(
                get    = '/session',
                PUTseal   = '/session/%s/seal',
                DELETE = '/session/%s'
            ),
            UPLOAD = list(
                uploadGenes = '/upload?session=%s&displayName=%s'
            ),
            JOB     = list(
                getAlgorithms = '/job/algorithms',
                getParameters = '/job/parameters/%s',
                startJob      = '/job/start?session=%s',
                setupJob      = '/job/setup/%s?session=%s',
                stopJob       = '/job/stop?session=%s'
            ),
            RESOURCE = list(
            )
        )
    ),
    private = list(
        GT_URL = 'https://genetrail2.bioinf.uni-sb.de/api',
        sessionId = NA
    )
)

## FLOW of client is normally:
# login
# start session
# upload file
# setup job
# start job
# once job is done, download the resources,
