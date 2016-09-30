#!/usr/bin/Rscript

suppressMessages( require( stringr ) )
suppressMessages( require( xlsx ) )

options( java.parameters = "-Xmx32g" )
.joptions( java.parameters = "-Xmx32g" )

.jinit()
.jaddClassPath( dir( '~/.local/java-libs', full.names = TRUE ) )
# Can prove it with .jclassPath()

# Java imports
    IndexedColors     <- 'org.apache.poi.ss.usermodel.IndexedColors'
    CellRangeAddress  <- 'org.apache.poi.ss.util.CellRangeAddress'
    CTAutoFilterImpl  <- 'org.openxmlformats.schemas.spreadsheetml.x2006.main.impl.CTAutoFilterImpl'
    STFilterOperator  <- 'org.openxmlformats.schemas.spreadsheetml.x2006.main.STFilterOperator'
    AreaReference     <- 'org.apache.poi.ss.util.AreaReference'
    CellReference     <- 'org.apache.poi.ss.util.CellReference'

# Excel table cell styles
    TABLE_ROWNAMES_STYLE <- list(
            fontName = 'Verdana', fontHeight = 12, fontIsBold = FALSE,
            alignment = Alignment( wrapText = TRUE, horizontal = "ALIGN_CENTER",
            vertical = "VERTICAL_CENTER" ) )
    TABLE_COLNAMES_STYLE <- list(
            fontName = 'Calibri', fontHeight = NULL, fontIsBold = TRUE,
            alignment = Alignment( wrapText = TRUE, horizontal = "ALIGN_CENTER",
            vertical = "VERTICAL_CENTER" ), dataFormat = DataFormat( '@' ) )
    TABLE_CELL_STYLE_TEXT <- list( list(
            fontName = 'Consolas', fontHeight = 12, fontIsBold = FALSE,
            alignment = Alignment( vertical = "VERTICAL_CENTER" ),
            dataFormat = DataFormat( '@' ) ) )
    TABLE_CELL_STYLE_NUMBER <- list( list(
            fontName = 'Consolas', fontHeight = 12, fontIsBold = FALSE,
            alignment = Alignment( vertical = "VERTICAL_CENTER" ),
            dataFormat = DataFormat( '0.00000' ) ) )
    TABLE_CELL_STYLE_NUMBER2 <- list( list(
            fontName = 'Consolas', fontHeight = 12, fontIsBold = FALSE,
            alignment = Alignment( vertical = "VERTICAL_CENTER" ),
            dataFormat = DataFormat( '0.00' ) ) )
    TABLE_CELL_STYLE_INTEGER <- list( list(
            fontName = 'Consolas', fontHeight = 12, fontIsBold = FALSE,
            alignment = Alignment( vertical = "VERTICAL_CENTER" ),
            dataFormat = DataFormat( '0' ) ) )
    TABLE_CELL_STYLE_EXP <- list( list(
            fontName = 'Consolas', fontHeight = 12, fontIsBold = FALSE,
            alignment = Alignment( vertical = "VERTICAL_CENTER" ),
            dataFormat = DataFormat( '0.0000E+00' ) ) )
    TABLE_CELL_STYLE_PERCENT <- list( list(
            fontName = 'Consolas', fontHeight = 12, fontIsBold = FALSE,
            alignment = Alignment( vertical = "VERTICAL_CENTER" ),
            dataFormat = DataFormat( '0.00%' ) ) )

    ENRICHED_COLOR <- J(IndexedColors)$ROSE$index
    DEPLETED_COLOR <- J(IndexedColors)$LIGHT_GREEN$index

##
# Create an object that can be stored as a sheet upon calling the outputToExcelFile function
#
createExcelSheetDataObject <- function(
    sheetName, dataFrame, cellStyles = c(), sheetTabColor = NULL,
    filter = NULL, columnWidths = 'auto' )
{
    if ( length( cellStyles ) != ncol( dataFrame ) )
    {
        if( length( cellStyles ) == 1 )
            # If there is only one style given, use that style for all columns
            cellStyles <- rep( cellStyles, ncol( dataFrame ) )
        else
        {
            logger( 'Number of cell styles does not match the number of data columns', level = logger.levels$ERROR )
            return()
        }
    }
    names( cellStyles ) <- 1:ncol( dataFrame )

    if ( length( columnWidths ) != ncol( dataFrame ) )
    {
        if( length( columnWidths ) == 1 )
            # If there is only one width given, use that width for all columns
            columnWidths <- rep( columnWidths, ncol( dataFrame ) )
        else
        {
            logger( 'Number of cell widths does not match the number of data columns', level = logger.levels$ERROR )
            return()
        }
    }

    return( list( list( 'sheetName' = sheetName, 'dataFrame' = dataFrame,
        'cellStyles' = cellStyles, 'sheetTabColor' = sheetTabColor,
        'filter' = filter, 'columnWidths' = columnWidths ) ) )
}

##
# Create an table filter object that is normally passed to the createExcelSheetDataObject function.
#  Useful if you wish to store more data in the excel file than would be helpful to actually show.
#
createTableFilter <- function(
    columnPattern = 'q_value',
    operator = J(STFilterOperator)$LESS_THAN_OR_EQUAL,
    val = 0.05 )
{
    return( list( 'columnPattern' = columnPattern, 'operator' = operator, 'val' = val ) )
}

##
# Write date to Excel file
#
outputToExcelFile <- function( allData = c(),
    fileName = 'secondaryAnalysisResults.xlsx' )
{
    logger( 'Creating Excel File', level = logger.levels$STAGE, append = '\n' )
    logger( fileName, level = logger.levels$FILE_PATH, append = '\n' )

    wb <- createWorkbook( type = 'xlsx' )

    setUpFilter <- function( tableFilter, data, sheet,
        applyFilterToRows = TRUE, columnPattern = 'q_value',
        operator = J(STFilterOperator)$LESS_THAN_OR_EQUAL, val = 0.05 )
    {
        fC <- grep( columnPattern, names( data ) ) - 1L
        filterColumn <- tableFilter$getFilterColumnArray( fC )
        customFilters <- filterColumn$addNewCustomFilters()
        customFilter <- customFilters$addNewCustomFilter()
        customFilter$setOperator( operator )
        customFilter$setVal( as.character( val ) )
        if ( applyFilterToRows )
            for ( row in getRows( sheet, -1 ) ) # Not the header row
            {
                cell <- row$getCell( fC )
                if ( getCellValue( cell ) >= val )
                    row$getCTRow()$setHidden( TRUE )
            }
    }

    sizeUpColumns <- function( sheet, widths )
    {
        for( i in 1:length( widths ) )
            if( !is.na( as.integer( widths[i] ) ) )
                setColumnWidth( sheet, i, as.integer( widths[i] ) )
            else autoSizeColumn( sheet,  i )
    }

    setUpTable <- function( data )
    {
        sheet <- createSheet( wb, sheetName = data$sheetName )
        if ( !is.null( data$sheetTabColor ) )
            sheet$setTabColor( data$sheetTabColor )

        areaReference <- new( J(AreaReference),
            new( J(CellReference), 0L, 0L ),
            new( J(CellReference), nrow(data$dataFrame) - 0L, ncol(data$dataFrame) - 1L ) )
        tab <- sheet$createTable()
        table <- tab$getCTTable()
        table_style <- table$addNewTableStyleInfo()
        table_style$setName( 'TableStyleLight1' )
        table_style$setShowRowStripes( TRUE )
        table$setRef( areaReference$formatAsString() )
        table$setDisplayName( data$sheetName )
        table$setName( paste0( data$sheetName, '_Test' ) )
        table$setId( .jlong( 1 ) )
        tableCols <- table$addNewTableColumns()
        tableCols$setCount( .jlong( ncol(data$dataFrame) ) )

        ctaf <- table$addNewAutoFilter()

        for ( col in 1:ncol( data$dataFrame ) )
        {
            tableCol <- tableCols$addNewTableColumn()
            tableCol$setName( names( data$dataFrame )[col] )
            tableCol$setId( .jlong( col ) )
            ctfc <- ctaf$addNewFilterColumn()
            ctfc$setColId( .jlong( col - 1 ) );
        }

        cellStyles <- data$cellStyles
        for ( csi in names( data$cellStyles ) )
        {
            cs <- data$cellStyles[[csi]]
            cellStyles[[csi]] <- CellStyle( wb,
                font = Font( wb, name = cs$fontName,
                    heightInPoints = cs$fontHeight, isBold = cs$fontIsBold ),
                alignment = cs$alignment,
                dataFormat = cs$dataFormat )
        }

        colNameCellStyle <- CellStyle( wb,
                font = Font( wb, name = TABLE_COLNAMES_STYLE$fontName,
                    isBold = TABLE_COLNAMES_STYLE$fontIsBold ),
                alignment = TABLE_COLNAMES_STYLE$alignment,
                dataFormat = TABLE_COLNAMES_STYLE$dataFormat )
        rowNameCellStyle <- CellStyle( wb,
                font = Font( wb, name = TABLE_ROWNAMES_STYLE$fontName,
                   isBold = TABLE_ROWNAMES_STYLE$fontIsBold ),
                alignment = TABLE_ROWNAMES_STYLE$alignment,
                dataFormat = TABLE_ROWNAMES_STYLE$dataFormat )
        addDataFrame( data$dataFrame,
            sheet = sheet, row.names = FALSE,
            startRow = 1, startColumn = 1,
            colnamesStyle = colNameCellStyle,
            rownamesStyle = rowNameCellStyle,
            colStyle = cellStyles )

        if ( !is.null( data$filter ) )
           setUpFilter( ctaf, data$dataFrame, sheet, columnPattern = data$filter$columnPattern,
               operator = data$filter$operator, val = data$filter$val )

        sizeUpColumns( sheet, data$columnWidths )

        return( table )
    }

    for ( data in allData )
        table <- setUpTable( data )

    saveWorkbook( wb, fileName )
}
