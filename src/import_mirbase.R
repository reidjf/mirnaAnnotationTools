#############################################################
## import_mirbase.R                                        ##
##  create annotation package for miRBase database         ##
##  using AnnotationDbi.                                   ##
##                                                         ##
## Author: James F. Reid <james.reid@ifom-ieo-campus.it>   ##
## Created: Mar 2010                                       ##
## Time-stamp: <2011-04-01 10:11:30 (jfr)>                 ##
##                                                         ##
## Notes: current working dir 'mirTools'                   ##
##        to update releases change:                       ##
##         dbVersion (package version of mirbase.db)       ##
##         mirbaseVersion                                  ##
##         mirbaseData                                     ##
##                                                         ##
#############################################################

#############################################################
## Table of Contents:
##  1. variable names and function definitions
##  2. download data from mirbase.org
##  3. prepare sqlite database (table definitions)
##  4. extract extra information not available in database
##  5. import data into database
##  6. build package by pointing to template directory
#############################################################

#############################################################
## variable names and paths
##
db                <- "mirbase"
package           <- sprintf("%s.db", db)
dbVersion         <- "0.5.0"
dbSchemaVersion   <- "2.1" ## AnnotationDbi dictates this
mirbaseVersion    <- "16.0"
mirbaseDate       <- "10 Sep 2010"
mirbaseDomain     <- "mirbase.org"
mirbaseCurrent    <- "pub/mirbase/CURRENT"
dbFTP             <- sprintf("ftp://%s/%s/", mirbaseDomain, mirbaseCurrent)
dbFTPdata         <- sprintf("%sdatabase_files/", dbFTP)
dbTitle           <- "miRBase: the microRNA database"
dbAnnObjTarget    <- "miRBase"
dbBiocViews       <- c("AnnotationData, FunctionalAnnotation")
dbAuthor          <- "James F. Reid <james.reid@ifom-ieo-campus.it>"
metaDataTables    <- c("metadata", "map_counts", "map_metadata")
dbOrganism        <- "Multiple"
dbSpecies         <- "Multiple"
dbManufacturer    <- "None"
dbManufacturerUrl <- "None"

## original data directory
dbRootPath      <- file.path("databases", db)
dbRootDataPath  <- file.path(dbRootPath, "data")

localdbFTP      <- file.path(dbRootPath, sub("ftp://", "", dbFTP))
localdbFTPdata  <- file.path(dbRootPath, sub("ftp://", "", dbFTPdata))

## template directory
templatePath <- file.path(dbRootPath, toupper(package))
dbSchemaFile <- file.path(dbRootPath, sprintf("%s_DB.sql", toupper(db)))
dbFile       <- file.path(tempdir(), sprintf("%s.sqlite", db))

## package directory
packagePath <- "packages"

## create directories if necc.
if (!file.exists(packagePath)) dir.create(packagePath)
if (!file.exists(dbRootDataPath)) dir.create(dbRootDataPath)
## note other ones MUST exist

downloadMirBase <- FALSE

#############################################################
## functions
##

## download all files from an ftp directory
downloadFtpPath <- function(url, destPath=tempdir()) {
    require("RCurl") || stop("Could not load package 'RCurl'.")

    downloadError <- function(e) stop("Could not download file.\n", e)

    ## get directory contents
    fileList <- tryCatch(getURL(url, ftp.use.epsv=FALSE),
                         error=downloadError)
    fileList <- strsplit(fileList, "\r*\n")[[1]]
    ## remove directories
    fileList <- fileList[sapply(fileList, function(i) {
        substr(strsplit(i, " ")[[1]][1], 1, 1)}) != "d"]
    ## names only
    fileList <- sapply(fileList, function(i) {
        tmp <- unlist(strsplit(i, " "))
        tmp[length(tmp)]})
    fileList <- sprintf("%s%s", url, fileList)

    ## download each file within directory
    for (ff in fileList) {
        filename <- basename(ff)
        tryCatch(download.file(url=ff,
                               destfile=file.path(destPath, filename),
                               mode="wb"), error=downloadError)
    }
}

## open table data file
readDataFile <- function(file) {
    return(read.delim(file=gzfile(file),
                      header=FALSE, stringsAsFactors=FALSE))
}

## write table data file
writeDataFile <- function(data, file) {
    con <- gzfile(file, open="w")
    write.table(data, file=con, sep="\t", row.names=FALSE, col.names=FALSE)
    close(con)
}

## insert data.frame 'data' into table 'dbTable' of database 'database'
sqlInsertTable <- function(database, dbTable, data) {
    require("RSQLite") || stop("Could not load package 'RSQLite'.")
    dbBeginTransaction(database)
    sqlString <- paste("INSERT INTO \'", dbTable, "\' VALUES(",
                 paste(":", colnames(data), sep="", collapse=","),
                       ");", sep="")
    dbGetPreparedQuery(database, sqlString, bind.data=data)
    dbCommit(database)
}

## count unique 'keys' in 'table' over connection 'con'
countDistinct <- function(con, table, key) {
    query <- sprintf("SELECT COUNT(DISTINCT %s) FROM %s;", key, table)
    return(sqliteQuickSQL(con, query)[1, 1])
}

#############################################################
## download data from miRBase ftp site
##
if (downloadMirBase) {
    cat("Downloading data...\n")
    dir.create(localdbFTPdata, recursive=TRUE)
    ## root files
    downloadFtpPath(url=dbFTP, destPath=localdbFTP)
    ## database files (MySQL data dumps)
    downloadFtpPath(url=dbFTPdata, destPath=localdbFTPdata)
}

#############################################################
## create extra tables
##

mirnaHairpinFile <- file.path(dbRootDataPath, "mirna_hairpin.txt.gz")
if (!file.exists(mirnaHairpinFile)) {
    cat("Creating extra 'mirna_hairpin' table.\n")
    ## add folding conformation of stem-loop sequence (not in db data)
    ## created with RNAfold program from the ViennaRNA suite.
    ## Hofacker IL, Stadler PF. Memory efficient folding algorithms for
    ## circular RNA secondary structures. Bioinformatics. 2006 May 15;
    ## 22(10):1172-6. PMID: 16452114
    miRNAstrFile <- file.path(localdbFTP, "miRNA.str.gz")
    con <- gzfile(miRNAstrFile)
    miRNAstr <- readLines(con)
    close(con)
    miRNAstr <- miRNAstr[which(miRNAstr != "")]
    ## scan each entry (six lines per sequence)
    strIndex <- grep(">", miRNAstr)
    if (!all(seq(1, length(miRNAstr), by=6) == strIndex)) {
        stop("Problem scanning ", miRNAstrFile)
    }
    tmp1 <- strsplit(miRNAstr[strIndex], " ")
    ## miRNA IDs
    rfIDs <- sub(">", "", unlist(lapply(tmp1, function(i) i[1])))
    ## minimum free energy
    rfMFE <- as.numeric(sub("\\(", "",
                            sub("\\)", "",
                                unlist(lapply(tmp1, function(i) i[2])))))
    ## mature/minor matches (this info is in mirna_mature table)
    ##tmp2 <- strsplit(miRNAstr[strIndex], "\\[")
    ##rfMMpos <- unlist(lapply(tmp2, function(i) {
    ##    paste(sub("\\]", "", sub(" ", "", i[2:length(i)])), collapse=", ")}))
    ## stem-loop folded sequence
    rfStemLoop <- sapply(strIndex, function(i) {
        paste(miRNAstr[(i+1):(i+5)], collapse="\n")})

    x <- data.frame('mirna_id' = rfIDs,
                    ##'mature_pos' = rfMMpos,
                    'hairpin' = rfStemLoop,
                    'mfe' = rfMFE,
                    check.names=FALSE, stringsAsFactors=FALSE)
    writeDataFile(x, mirnaHairpinFile)
}

mirnaClusterFile <- file.path(dbRootDataPath, "mirna_cluster.txt.gz")
if (!file.exists(mirnaClusterFile)) {
    cat("Creating extra 'mirna_cluster' table.\n")
    ## search for 'clustered' mirna in each species
    ## ie. other mirnas <10kb from any given mirna
    clusterWindow <- 10000

    ## read-in genomic coordinates
    coords <- readDataFile(file.path(localdbFTPdata,
                                     "mirna_chromosome_build.txt.gz"))
    colnames(coords) <- c('_id', 'xsome', 'contig_start',
                          'contig_end', 'strand')
    ## read-in mirna ids
    mirna <- readDataFile(file.path(localdbFTPdata, "mirna.txt.gz"))
    colnames(mirna) <- c('_id', 'mirna_acc', 'mirna_id', 'description',
                         'sequence', 'comment', 'auto_species')

    mirnaCoords <- cbind(mirna[match(coords[, '_id'], mirna[, '_id']),
                               c(1, 3, 7)], coords[, 2:5])

    ## initialize mirna_cluster table and cluster id
    mirnaCluster <- data.frame()
    ## iterate over each organism (48 with genomic coordinates - v.14)
    for (orgID in unique(mirnaCoords$auto_species)) {
        cat("\n", orgID, "\n")
        ## list of miRNAs belonging to this species
        mirnaID <- mirnaCoords[mirnaCoords$auto_species == orgID, 'mirna_id']
        ## compute cluster window around each member

        ## extract chromosome name
        mirnaChr <- mirnaCoords[(mirnaCoords[, 'mirna_id'] %in% mirnaID),
                                c('mirna_id', 'xsome')]
        ## remove multiple mappings (mirbase.org does NOT do this...)
        mult <- unique(mirnaChr[duplicated(mirnaChr[, 1]), 1])
        mirnaChr <- mirnaChr[!(mirnaChr[, 'mirna_id'] %in% mult), ]

        ## skip loop if nothing is left
        if (nrow(mirnaChr) == 0) next

        ## append start and end of mirna
        mcMatch <- (mirnaCoords[, 'mirna_id'] %in% mirnaChr[, 'mirna_id'])
        mirnaChr <- cbind(mirnaChr,
                          mirnaCoords[mcMatch, 'contig_start'],
                          mirnaCoords[mcMatch, 'contig_end'])
        ## and start (s) and end (e) of cluster window
        ## (note: strand is not taken into account)
        mirnaChr <- cbind(mirnaChr,
                          c(mirnaChr[, 3] - clusterWindow),
                          c(mirnaChr[, 4] + clusterWindow))
        colnames(mirnaChr) <- c('mirna_id', 'chromosome',
                                'cs', 'ce', 's', 'e')
        ## order by chromosome then by start
        mirnaChr <- with(mirnaChr, mirnaChr[order(mirnaChr[, 'chromosome'],
                                                  mirnaChr[, 's']), ])
        ## scan through each chromosome
        for (chr in unique(mirnaChr$chromosome)) {
            cat("\t", chr, ".")
            tmp <- mirnaChr[mirnaChr$chromosome == chr, ]
            if (nrow(tmp) == 1) next
            ## overlap matches
            tmpM <- sapply(tmp$cs, function(mir) {
                (mir > tmp$s & mir < tmp$e)})
            rownames(tmpM) <- colnames(tmpM) <- tmp[, 'mirna_id']
            ## collect results
            clMember <- lapply(1:nrow(tmpM), function(i) {
                rownames(tmpM)[tmpM[, i]]})
            clCluster <- lapply(1:nrow(tmpM), function(i) {
                rep(i, length(clMember[[i]]))})
            clId <- lapply(1:nrow(tmpM), function(i) {
                rep(rownames(tmpM)[i], length(clMember[[i]]))})
            mirnaCluster <- rbind(mirnaCluster,
                                  cbind(unlist(clId),
                                        unlist(clMember),
                                        unlist(clCluster)))
        }
    }
    colnames(mirnaCluster) <- c('mirna_id', 'member', 'cluster')
    writeDataFile(mirnaCluster, mirnaClusterFile)
}

#############################################################
## modify original tables
##
cat("Modifying original data tables\n")

dataFiles <- c("literature_references.txt.gz",
               "mirna_2_prefam.txt.gz",
               "mirna_chromosome_build.txt.gz",
               "mirna_context.txt.gz",
               "mirna_database_links.txt.gz",
               "mirna_literature_references.txt.gz",
               "mirna_mature.txt.gz",
               "mirna_prefam.txt.gz",
               "mirna_pre_mature.txt.gz",
               "mirna_species.txt.gz",
               "mirna.txt.gz")

for (ff in dataFiles) {
    if (ff == "literature_references.txt.gz") {
        cat("Replacing '\\N' with NA in 'literature_references' table...\n")
        x <- readDataFile(file.path(localdbFTPdata, ff))
        colnames(x) <- c("auto_lit", "medline", "title", "author", "journal")
        ## substitute "\\N" with NA
        x[x[, 'medline'] == "\\N", 'medline'] <- NA
        writeDataFile(x, file.path(dbRootDataPath, ff))
    }
    if (ff == "mirna_2_prefam.txt.gz") {
        cat("Copying 'mirna_2_prefam' table...\n")
        file.copy(file.path(localdbFTPdata, ff), file.path(dbRootDataPath, ff))
    }
    if (ff == "mirna_chromosome_build.txt.gz") {
        cat("Appending strand information to genomic coordinates of",
            "'mirna_chromosome_build' table...\n")
        x <- readDataFile(file.path(localdbFTPdata, ff))
        colnames(x) <- c("_id", "xsome", "contig_start",
                         "contig_end", "strand")
        ## incorporate strand (+/-) to coordinates to conform with CHR,
        ## CHRLOC and CHRLOCEND of other BioC annotation packages
        strand <- ifelse(x[, 'strand'] == "-", -1, 1)
        x[, 'contig_start'] <- strand * x[, 'contig_start']
        x[, 'contig_end'] <- strand * x[, 'contig_end']
        writeDataFile(x, file.path(dbRootDataPath, ff))
    }
    if (ff == "mirna_context.txt.gz") {
        cat("Copying 'mirna_context' table...\n")
        file.copy(file.path(localdbFTPdata, ff), file.path(dbRootDataPath, ff))
    }
    if (ff == "mirna_database_links.txt.gz") {
        cat("Dropping empty columns in 'mirna_database_links' table...\n")
        x <- readDataFile(file.path(localdbFTPdata, ff))
        colnames(x) <- c("_id", "db_id", "comment", "db_link",
                         "db_secondary", "other_params")
        ## remove target links (targetscan, pictar, etc.)
        ## "TARGETS:MIRTE"       "TARGETS:PICTAR-FLY"  "TARGETS:PICTAR-VERT"
        targetZ <- grep("TARGETS", levels(factor(x[, 'db_id'])))
        ## and "MIRTE"
        targetZ2 <- unique(c(targetZ,
                             grep("MIRTE", levels(factor(x[, 'db_id'])))))
        x <- x[-which(x[, 'db_id'] %in%
                      levels(factor(x[, 'db_id']))[targetZ2]), ]
        ## substitute "\\N" with NA
        x[x[, 'db_secondary'] == "\\N", 'db_secondary'] <- NA
        ## remove empty columns 'comment' and 'other-params'
        x <- x[,  c("_id", "db_id", "db_link", "db_secondary")]
        writeDataFile(x, file.path(dbRootDataPath, ff))
    }
    if (ff == "mirna_literature_references.txt.gz") {
        cat("Dropping empty columns in 'mirna_literature_references'",
            "table...\n")
        x <- readDataFile(file.path(localdbFTPdata, ff))
        colnames(x) <- c("_id", "auto_lit", "comment", "order_added")
        ## remove empty column 'comment'
        x <- x[,  c("_id", "auto_lit", "order_added")]
        writeDataFile(x, file.path(dbRootDataPath, ff))
    }
    if (ff == "mirna_mature.txt.gz") {
        cat("Replacing '\\N' with NA in 'mirna_mature' table...\n")
        x <- readDataFile(file.path(localdbFTPdata, ff))
        colnames(x) <- c("auto_mature", "mature_name", "mature_acc",
                         "mature_from", "mature_to", "evidence", "experiment",
                         "similarity")
        ## substitute "\\N" with NA
        x[x[, 'experiment'] == "\\N", 'experiment'] <- NA
        writeDataFile(x, file.path(dbRootDataPath, ff))
    }
    if (ff == "mirna_prefam.txt.gz") {
        cat("Dropping empty columns in 'mirna_prefam' table...\n")
        x <- readDataFile(file.path(localdbFTPdata, ff))
        colnames(x) <- c("auto_prefam", "prefam_acc", "prefam_id",
                         "description")
        ## remove empty column 'description'
        x <- x[,  c("auto_prefam", "prefam_acc", "prefam_id")]
        writeDataFile(x, file.path(dbRootDataPath, ff))
    }
    if (ff == "mirna_pre_mature.txt.gz") {
        cat("Copying 'mirna_pre_mature' table...\n")
        file.copy(file.path(localdbFTPdata, ff), file.path(dbRootDataPath, ff))
    }
    ## species file handled together with 'mirna'
    if (ff == "mirna.txt.gz") {
        cat("Appending 'species' to mirna table...\n")
        mirna <- readDataFile(file.path(localdbFTPdata, ff))
        colnames(mirna) <- c("_id", "mirna_acc", "mirna_id", "description",
                             "sequence", "comment", "auto_species")
        mirnaSpecies <- readDataFile(file.path(localdbFTPdata,
                                               "mirna_species.txt.gz"))
        colnames(mirnaSpecies) <- c("auto_species", "organism", "division",
                                    "name", "taxonomy", "genome_assembly",
                                    "ensembl_db")
        mirnaOrg <- mirnaSpecies[match(mirna$auto_species,
                                       mirnaSpecies$auto_species), 'organism']
        mirna <- cbind(mirna[, 1:(ncol(mirna) - 1)], organism=mirnaOrg)

        ## replace empty strings with NA in 'comment' column
        ##mirna[mirna$comment == "", 'comment'] <- NA
        writeDataFile(mirna, file.path(dbRootDataPath, ff))
        writeDataFile(mirnaSpecies, file.path(dbRootDataPath,
                                              "mirna_species.txt.gz"))
    }
}
## we don't import mirna_target_links.txt.gz and mirna_target_url.txt.gz


#############################################################
## create sqlite db
##
cat("Create SQL database...\n")

## initialize database
unlink(dbFile)
require("RSQLite") || stop("Could not load package 'RSQLite'.")
drv <- dbDriver("SQLite")
con <- dbConnect(drv, dbname=dbFile)

## read db schema
dbSchema <- readLines(dbSchemaFile)
dbSchema <- strsplit(paste(dbSchema, collapse="\n"), ";")[[1]]
## create tables
dbTablesSQL <- sapply(dbSchema[grep("CREATE TABLE", dbSchema)],
                      function(x) sqliteQuickSQL(con, x))
## extract table names
dbTables <- unlist(lapply(names(dbTablesSQL), function(x) {
    strsplit(strsplit(x, "CREATE TABLE ")[[1]][2], " \\(\\n")[[1]][1]}))
## remove meta-data tables
dbTables <- dbTables[!(dbTables %in% metaDataTables)]

## store number of elements in each table
mbT <- vector()

## loop over each table
for (tt in dbTables) {
    x <- read.delim(file=gzfile(file.path(dbRootDataPath,
                    sprintf("%s.txt.gz", tt))),
                    header=FALSE, stringsAsFactors=FALSE)
    colnames(x) <- sqliteTableFields(con, tt)
    ## insert data
    sqlInsertTable(con, tt, x)
    mbT[tt] <- nrow(x)
}

## append the metadata info

## metadata table
metaData <- rbind(c("DBSCHEMA", sprintf("%s_DB", toupper(db))),
                  c("ORGANISM", dbOrganism),
                  c("SPECIES", dbSpecies),
                  c("DBSCHEMAVERSION", dbSchemaVersion),
                  c(sprintf("%sSOURCENAME", toupper(db)), dbAnnObjTarget),
                  c(sprintf("%sSOURCEURL", toupper(db)), dbFTP),
                  c(sprintf("%sSOURCEDATE", toupper(db)), mirbaseDate),
                  ## JFR non-conventional
                  c(sprintf("%sSOURCEVERSION", toupper(db)), mirbaseVersion),
                  c(sprintf("%sVERSION", toupper(db)), dbVersion))
q <- paste("INSERT INTO 'metadata' VALUES('",
           metaData[, 1], "','", metaData[, 2], "');", sep="")
tmp <- sapply(q, function(x) sqliteQuickSQL(con, x))

## map_counts table
mapCounts <- rbind(c("CHR",
                     countDistinct(con, "mirna_chromosome_build", "_id")),
                   c("CHRLOC",
                     countDistinct(con, "mirna_chromosome_build", "_id")),
                   c("CHRLOCEND",
                     countDistinct(con, "mirna_chromosome_build", "_id")),
                   c("CLUSTER",
                     countDistinct(con, "mirna_cluster", "mirna_id")),
                   c("COMMENT",
                     countDistinct(con, "mirna", "mirna_id")),
                   c("CONTEXT",
                     countDistinct(con, "mirna_context", "_id")),
                   c("DESCRIPTION",
                     countDistinct(con, "mirna", "mirna_id")),
                   c("FAMILY",
                     countDistinct(con, "mirna_2_prefam", "_id")),
                   c("HAIRPIN",
                     countDistinct(con, "mirna_hairpin", "mirna_id")),
                   c("ID2ACC",
                     countDistinct(con, "mirna", "_id")),
                   c("ID2SPECIES",
                     countDistinct(con, "mirna", "_id")),
                   c("LINKS",
                     countDistinct(con , "mirna_database_links", "_id")),
                   c("MATURE",
                     countDistinct(con, "mirna_pre_mature", "_id")),
                   c("MFE",
                     countDistinct(con, "mirna_hairpin", "mirna_id")),
                   c("PMID",
                     countDistinct(con, "mirna_literature_references", "_id")),
                   c("SEQUENCE",
                     countDistinct(con, "mirna", "_id")),
                   c("SPECIES",
                     countDistinct(con, "mirna_species", "organism")))
q <- paste("INSERT INTO 'map_counts' VALUES('",
           mapCounts[, 1], "',", mapCounts[, 2], ");", sep="")
tmp <- sapply(q, function(x) sqliteQuickSQL(con, x))

## map_metadata
mapMetadata <- data.frame(mapCounts[, 1],
                          sprintf("%s (Version: %s)", dbAnnObjTarget,
                                  mirbaseVersion), dbFTP, mirbaseDate)
q <-  paste("INSERT INTO 'map_metadata' VALUES('",
            mapMetadata[, 1], "', '", mapMetadata[, 2], "', '",
            mapMetadata[, 3], "', '", mapMetadata[, 4], "');", sep="")
tmp <- sapply(q, function(x) sqliteQuickSQL(con, x))

## sanity checks
for (i in 1:length(mbT)) {
    if (dbGetQuery(con,paste("SELECT COUNT(*) FROM ",
                             names(mbT)[i], sep="")) != mbT[i]) {
        stop("Element count mismatch.")
    }
}

## clean and disconnect
dbGetQuery(con, "VACUUM")
dbDisconnect(con)


#############################################################
## create package
##
cat("Creating package...\n")

## small modification to AnnotationDbi:::makeAnnDbPkg to handle
## template directory not under AnnotationDbi system folder
##library("AnnotationDbi")
source("src/makeAnnDbPkg_dev.R")

## build new package
seed <- new("AnnDbPkgSeed",
            Package = package,
            Title = dbTitle,
            Version = dbVersion,
            License="file LICENSE",
            Author = dbAuthor,
            Maintainer = dbAuthor,
            PkgTemplate = templatePath,
            DBschema = sprintf("%s_DB", toupper(db)),
            AnnObjPrefix = db,
            AnnObjTarget = dbAnnObjTarget,
            organism = dbOrganism,
            species = dbSpecies,
            manufacturer = dbManufacturer,
            ##chipName = as.character(NA)
            manufacturerUrl = dbManufacturerUrl,
            biocViews = dbBiocViews)

unlink(file.path(packagePath, package), recursive=TRUE)
build <- try(makeAnnDbPkg(seed, dbFile, dest_dir=packagePath, no.man=FALSE))

if (!inherits(build, 'try-error')) { ## install new package and load it!
    cat("Installing package...\n")
    install.packages(file.path(packagePath, package), repos=NULL, type="src")
    library(package, character.only=TRUE)
} else {
    unlink(file.path(packagePath, package), recursive=TRUE)
    unlink(dbFile)
}

#############################################################
## HISTORY:
## 2010-03-10
## o finished creation of mirbase.db annotation package
## 2010-03-24
## o added extensive man pages and methods for complex
## objects (CONTEXT, LINKS, MATURE and PMID)
## 2010-04-25
## o updated to mirbase v. 15
## 2010-09-30
## o updated to mirbae v. 16
## 2011-04-01
## o new bioc release, minor layout in BASE.Rd
#############################################################




