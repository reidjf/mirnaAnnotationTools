#############################################################
## import_targetscan.Mm.eg.R                               ##
##  create annotation package for targetscan (Mouse)       ##
##  database using AnnotationDbi.                          ##
##                                                         ##
## Author: James F. Reid <james.reid@ifom-ieo-campus.it>   ##
## Created: Mar 2010                                       ##
## Time-stamp: <2011-04-01 10:42:22 (jfr)>                 ##
##                                                         ##
## Notes: current working dir 'mirTools'                   ##
#############################################################

#############################################################
## Table of Contents:
##  1. variable names and function definitions
##  2. download data from targetscan.org
##  3. prepare sqlite database (table definitions)
##  4. extract extra information not available in database
##  5. import data into database
##  6. build package by pointing to template directory
#############################################################

#############################################################
## variable names and paths
##
db                 <- "targetscan.Mm.eg"
package            <- sprintf("%s.db", db)
dbVersion          <- "0.3.0"
dbSchemaVersion    <- "2.1" ## AnnotationDbi dictates this
targetscanVersion  <- "5.1"
targetscanDate     <- "April 2009"
targetscanDomain   <- "www.targetscan.org"
targetscanCurrent1 <- "mmu_50/mmu_50_data_download"
targetscanCurrent2 <- "cgi-bin/targetscan/data_download.cgi?db=vert_50"
dbFTP              <- sprintf("http://%s/%s/", targetscanDomain,
                              targetscanCurrent1)
dbTitle            <- "TargetScan miRNA target predictions for mouse"
dbAnnObjTarget     <- "TargetScan (Mouse)"
dbBiocViews        <- c("AnnotationData, FunctionalAnnotation")
dbAuthor           <- "Gabor Csardi <Gabor.Csardi@unil.ch>"
dbMaintainer       <- "James F. Reid <james.reid@ifom-ieo-campus.it>"
metaDataTables     <- c("metadata", "map_counts", "map_metadata")
dbOrganism         <- "Mus musculus"
dbSpecies          <- "Mouse"
dbManufacturer     <- "None"
dbManufacturerUrl  <- "None"

## original data directory
dbRootPath <- file.path("databases", db)
localdbFTP <- file.path(dbRootPath, sub("http://", "", dbFTP))

## template directory
templatePath  <- file.path(dbRootPath, toupper(package))
dbSchemaFile  <- file.path(dbRootPath, sprintf("%s_DB.sql", toupper(db)))
dbFile        <- file.path(tempdir(), sprintf("%s.sqlite", db))

## package directory
packagePath <- "packages"

## create directories if necc.
if (!file.exists(packagePath)) dir.create(packagePath)
if (!file.exists(dbRootPath)) dir.create(dbRootPath)

## note other ones MUST exist

#############################################################
## functions
##

downloadError <- function(e) stop("Could not download file.\n", e)

downloadArchive <- function(URL, file, destPath=tempdir()) {

    zipFile <- sprintf("%s.zip", file)
    URL <- sprintf("%s%s", URL, zipFile)
    tmpZipFile <- file.path(tempdir(), zipFile)
    tmpFile <- file.path(tempdir(), file)

    tryCatch(download.file(url=URL, destfile=tmpZipFile), mode="wb",
             error=downloadError)
    zip.file.extract(tmpFile, zipFile)
    file.copy(tmpFile, file.path(destPath, file))
}

#############################################################
## download data from targetscan and ncbi sites
##
targetfile <- file.path(localdbFTP, "Predicted_Targets_Info.txt")
familyfile <- file.path(localdbFTP, "miR_Family_Info.txt")
tgFiles <- c(basename(targetfile), basename(familyfile))
taxfile <- file.path(localdbFTP, "names.dmp")

if (!file.exists(dirname(localdbFTP))) {
    cat("Downloading data...\n")
    dir.create(localdbFTP, recursive=TRUE)

    for (ff in tgFiles) {
        cat("Downloading", ff, "...\n")
        downloadArchive(dbFTP, ff, localdbFTP)
    }
    taxfile.zip <- file.path(localdbFTP, "taxdmp.zip")
    if (!file.exists(taxfile)) {
        cat("Downloading", taxfile, "...\n")
        data.url <- "ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdmp.zip"
        tryCatch(download.file(data.url, destfile=taxfile.zip),
                          error=downloadError)
        taxfile.tmp <- zip.file.extract(taxfile,
                                        basename(taxfile.zip))
        file.copy(taxfile.tmp, taxfile)
        file.remove(taxfile.zip)
    }
}

dbFTP <- sprintf("http://%s/%s/", targetscanDomain, targetscanCurrent2)

#############################################################
## read in data
##
cat("Reading in data...\n")

family <- read.delim(familyfile)
targets <- read.delim(targetfile)

tax <- read.delim(taxfile, header=FALSE, quote="")
tax <- tax[, -c(2, 4, 6, 8)]
sci <- which(tax[, 4]=="scientific name")
tax <- tax[sci,]

species <- unique(family$Species.ID)
names <- tax[, 2][ match(species, tax[, 1]) ]
species <- data.frame(id=species, name=names)

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


###################################################
### put the data into the tables
###################################################
## Append data
## Species
dbBeginTransaction(con)
dbGetPreparedQuery(con, 'INSERT INTO "species" VALUES(:id, :name);',
                   bind.data=species)
dbCommit(con)

## miRNA families
family$miR.family <- sub("^miR-141/200$", "miR-141/200a",
                         family$miR.family, fixed=FALSE)
family$miR.family <- sub("^miR-200bc/420$", "miR-200bc/429",
                         family$miR.family, fixed=FALSE)
fam <- unique(as.character(family$miR.family))
fam <- cbind(id=seq_along(fam), fam)
dbBeginTransaction(con)
dbGetPreparedQuery(con,
                   "INSERT INTO 'mirna_family' VALUES(:id, :fam);",
                   bind.data=as.data.frame(fam))
dbCommit(con)

## mirbase table
mirbase <- family[,c("MiRBase.ID", "MiRBase.Accession",
                     "miR.family", "Seed.m8", "Species.ID",
                     "Mature.sequence",
                     "Family.Conservation.")]
mirbase$miR.family <- fam[,1][ match(family$miR.family, fam[,2]) ]
colnames(mirbase) <- letters[1:7]
dbBeginTransaction(con)
dbGetPreparedQuery(con, bind.data=mirbase,
                   'INSERT INTO "mirbase" VALUES(:a,:b,:c,:d,:e,:f,:g)')
dbCommit(con)

## keep entries for mouse only
targets2 <- targets
targets2 <- targets2[ targets2$Species.ID == 10090, ]

## genes
gs <- unique(targets2$Gene.ID)
gs <- cbind(seq_along(gs), gs)
dbBeginTransaction(con)
dbGetPreparedQuery(con, bind.data=data.frame(a=gs[,1],
                         b=as.integer(gs[,2])),
                   "INSERT INTO 'genes' VALUES(:a,:b);")
dbCommit(con)

## seed_match
sm <- sort(unique(as.character(targets$Seed.match)))
sm <- cbind(seq_along(sm), sm)
dbBeginTransaction(con)
dbGetPreparedQuery(con, bind.data=data.frame(a=sm[,1], b=sm[,2]),
                   "INSERT INTO 'seed_match' VALUES(:a,:b);")
dbCommit(con)

## targets, human only :(
targets2$miR.Family <- fam[,1][ match(targets2$miR.Family, fam[,2]) ]
targets2$Gene.ID <- gs[,1][ match(targets2$Gene.ID, gs[,2]) ]
targets2$Seed.match <- sm[,1][ match(targets2$Seed.match, sm[,2]) ]
colnames(targets2) <- sub(".", "_", colnames(targets2), fixed=TRUE)
dbBeginTransaction(con)
dbGetPreparedQuery(con, bind.data=targets2,
                   paste(sep="",
                         "INSERT INTO targets VALUES(:miR_Family,",
                         ":Gene_ID, :Species_ID,",
                         ":UTR_start, :UTR_end,",
                         ":MSA_start, :MSA_end,",
                         ":Seed_match, :PCT);"))
dbCommit(con)


## append the metadata info

## metadata table
metaData <- rbind(c("DBSCHEMA", sprintf("%s_DB", toupper(db))),
                  c("ORGANISM", dbOrganism),
                  c("SPECIES", dbSpecies),
                  c("DBSCHEMAVERSION", dbSchemaVersion),
                  c(sprintf("%sSOURCENAME", toupper(db)), dbAnnObjTarget),
                  c(sprintf("%sSOURCEURL", toupper(db)), dbFTP),
                  c(sprintf("%sSOURCEDATE", toupper(db)), targetscanDate),
                  ## JFR non-conventional
                  c(sprintf("%sSOURCEVERSION", toupper(db)), targetscanVersion),
                  c(sprintf("%sVERSION", toupper(db)), dbVersion))
q <- paste("INSERT INTO 'metadata' VALUES('",
           metaData[, 1], "','", metaData[, 2], "');", sep="")
tmp <- sapply(q, function(x) sqliteQuickSQL(con, x))

## map_counts table
mapCounts <- rbind(c("FAMILY2MIRBASE", nrow(fam)),
                   c("MIRBASE2FAMILY", nrow(mirbase)),
                   c("MIRNA", nrow(mirbase)),
                   c("TARGETS", nrow(gs)),
                   c("TARGETSFULL", nrow(gs))
                   )

q <- paste("INSERT INTO 'map_counts' VALUES('",
           mapCounts[, 1], "',", mapCounts[, 2], ");", sep="")
tmp <- sapply(q, function(x) sqliteQuickSQL(con, x))

## map_metadata
mapMetadata <- data.frame(mapCounts[, 1],
                          sprintf("%s (Version: %s)", dbAnnObjTarget,
                                  targetscanVersion), dbFTP, targetscanDate)
q <-  paste("INSERT INTO 'map_metadata' VALUES('",
            mapMetadata[, 1], "', '", mapMetadata[, 2], "', '",
            mapMetadata[, 3], "', '", mapMetadata[, 4], "');", sep="")
tmp <- sapply(q, function(x) sqliteQuickSQL(con, x))

## sanity checks
if (dbGetQuery(con, "SELECT COUNT(*) FROM species") != nrow(species)) {
  stop("FOOBAR")
}
if (dbGetQuery(con, "SELECT COUNT(*) FROM mirna_family") != nrow(fam)) {
  stop("FOOBAR")
}
if (dbGetQuery(con, "SELECT COUNT(*) FROM mirbase") != nrow(mirbase)) {
  stop("FOOBAR")
}
if (dbGetQuery(con, "SELECT COUNT(*) FROM genes") != nrow(gs)) {
  stop("FOOBAR")
}
if (dbGetQuery(con, "SELECT COUNT(*) FROM seed_match") != nrow(sm)) {
  stop("FOOBAR")
}
if (dbGetQuery(con, "SELECT COUNT(*) FROM targets") != nrow(targets2)) {
  stop("FOOBAR")
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
            Maintainer = dbMaintainer,
            PkgTemplate = templatePath,
            DBschema = sprintf("%s_DB", toupper(db)),
            AnnObjPrefix = db,
            AnnObjTarget = dbAnnObjTarget,
            organism = dbOrganism,
            species = dbSpecies,
            manufacturer = dbManufacturer,
            ##chipName = as.character(NA)
            manufacturerUrl = dbManufacturerUrl,
            biocViews = dbBiocViews
            )

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
## 2010-03-24
## o modifying original script written by Gabor to
##   conform to mirbase creation
## 2010-09-30
## o new release of BioC (no data change)
## 2011-04-01
## o new release of BioC (no data change)
## o added extra info viz type of targets available
## o added ref to mirbase.db, changed mirbase url
#############################################################





