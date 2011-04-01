datacache <- new.env(hash=TRUE, parent=emptyenv())

## append mirbase version and date
@ANNOBJPREFIX@ <- function() {
    showQCData("@ANNOBJPREFIX@", datacache)
    cat(paste("miRBase version: ",
              dbmeta(datacache, 'MIRBASESOURCEVERSION'), "\n", sep=""))
    cat(paste("miRBase date: ",
              dbmeta(datacache, 'MIRBASESOURCEDATE'), "\n", sep=""))
}
@ANNOBJPREFIX@_dbconn <- function() dbconn(datacache)
@ANNOBJPREFIX@_dbfile <- function() dbfile(datacache)
## FIXME(jfr) not working
##dbschema(mirbase_dbconn(), mirbase_dbfile())
##
##Warning messages:
##1: In file(con, "r") :
##  file("") only supports open = "w+" and open = "w+b": using the former
##2: In max(ii) : no non-missing arguments to max; returning -Inf
##@ANNOBJPREFIX@_dbschema <- function(file="", show.indices=FALSE) {
##    dbschema(datacache, file=file, show.indices=show.indices)
##}
## alternative by reading sql file directly
@ANNOBJPREFIX@_dbschema <- function(file="", show.indices=FALSE) {
    x <- dbconn(datacache)
    q <- paste("SELECT sql FROM (SELECT * FROM sqlite_master UNION ALL",
               "SELECT * FROM sqlite_temp_master) WHERE type!='meta'",
               "ORDER BY tbl_name, type DESC, name", sep=" ")
    schema <- sqliteQuickSQL(x, q)
    cat(schema[!is.na(schema), ], sep=";\n\n")
}

@ANNOBJPREFIX@_dbInfo <- function() dbInfo(datacache)


.onLoad <- function(libname, pkgname) {
    require("methods", quietly=TRUE)

    ## custom Bimap classes
    ## 'complex' Bimaps are 'CONTEXT', 'LINKS', 'MATURE' and 'PMID'
    ## each Bimap named 'NAME' is called 'mirbaseNameBimap' and has an
    ## associated class for representing the output 'mirnaName'

    ## CONTEXT
    setClass("mirbaseContextBimap", contains="AnnDbBimap")
    setClass("mirnaContext",
             representation(contextTranscriptID="character",
                            contextOverlapSense="character",
                            contextOverlapType="character",
                            contextNumber="numeric",
                            contextTranscriptSource="character",
                            contextTranscriptName="character"))

    setMethod("as.list", "mirbaseContextBimap", function(x, ...) {
        y <- AnnotationDbi:::flatten(x, fromKeys.only=TRUE)
        makeMBcontextNode <- function(transcript_id, overlap_sense,
                                      overlap_type, number, transcript_source,
                                      transcript_name, ...) {
            new("mirnaContext",
                contextTranscriptID=transcript_id,
                contextOverlapSense=overlap_sense,
                contextOverlapType=overlap_type, contextNumber=number,
                contextTranscriptSource=transcript_source,
                contextTranscriptName=transcript_name, ...)
        }
        AnnotationDbi:::.toListOfLists(y, mode=1, makeMBcontextNode)
    })
    setGeneric("contextTranscriptID",
               function(object) standardGeneric("contextTranscriptID"))
    setMethod("contextTranscriptID", "mirnaContext",
              function(object) object@contextTranscriptID)
    setGeneric("contextOverlapSense",
               function(object) standardGeneric("contextOverlapSense"))
    setMethod("contextOverlapSense", "mirnaContext",
              function(object) object@contextOverlapSense)
    setGeneric("contextOverlapType",
               function(object) standardGeneric("contextOverlapType"))
    setMethod("contextOverlapType", "mirnaContext",
              function(object) object@contextOverlapType)
    setGeneric("contextNumber",
               function(object) standardGeneric("contextNumber"))
    setMethod("contextNumber", "mirnaContext",
              function(object) object@contextNumber)
    setGeneric("contextTranscriptSource",
               function(object) standardGeneric("contextTranscriptSource"))
    setMethod("contextTranscriptSource", "mirnaContext",
              function(object) object@contextTranscriptSource)
    setGeneric("contextTranscriptName",
               function(object) standardGeneric("contextTranscriptName"))
    setMethod("contextTranscriptName", "mirnaContext",
              function(object) object@contextTranscriptName)

    setMethod("show", "mirnaContext", function(object) {
        s <- character()
        strand <- object@contextOverlapSense
        plusMatch <- strand %in% "+"
        minusMatch <- strand %in% "-"
        printContext <- function(object, index) {
            pc <- sprintf("  %s; %s; %s (%s) [%s]",
                          object@contextTranscriptID[index],
                          object@contextTranscriptName[index],
                          object@contextOverlapType[index],
                          object@contextNumber[index],
                          object@contextTranscriptSource[index])
            return(sub("; ; ", "; ", sub("\\[\\]", "", pc)))
        }

        if (sum(plusMatch) > 0) {
            s <- c(s, cat("[+]\n"))
            s <- c(s, printContext(object, plusMatch))
            cat(s, sep="\n")
        }
        s <- character()
        if (sum(minusMatch) > 0) {
            s <- c(s, cat("[-]\n"))
            s <- c(s, printContext(object, minusMatch))
            cat(s, sep="\n")
        }
    })

    ## LINKS
    setClass("mirbaseLinksBimap", contains="AnnDbBimap")
    setClass("mirnaLinks",
             representation(linksDbLink = "character",
                            linksDbId = "character",
                            linksDbSecondary = "character"))

    setMethod("as.list", "mirbaseLinksBimap", function(x, ...) {
        y <- AnnotationDbi:::flatten(x, fromKeys.only=TRUE)
        makeMBlinksNode <- function(db_link, db_id, db_secondary, ...) {
            new("mirnaLinks",
                linksDbLink=db_link, linksDbId = db_id,
                linksDbSecondary = db_secondary, ...)
        }
        AnnotationDbi:::.toListOfLists(y, mode=1, makeMBlinksNode)
    })
    setGeneric("linksDbLink", function(object) standardGeneric("linksDbLink"))
    setMethod("linksDbLink", "mirnaLinks", function(object) object@linksDbLink)
    setGeneric("linksDbId", function(object) standardGeneric("linksDbId"))
    setMethod("linksDbId", "mirnaLinks", function(object) object@linksDbId)
    setGeneric("linksDbSecondary",
               function(object) standardGeneric("linksDbSecondary"))
    setMethod("linksDbSecondary", "mirnaLinks",
              function(object) object@linksDbSecondary)

    setMethod("show", "mirnaLinks", function(object) {
        s <- character()
        dbNames <- unique(object@linksDbId)
        IDnames <- sub("\\(NA\\)", "",
                       paste("(", object@linksDbSecondary, ")", sep=""))
        IDlinks <- object@linksDbLink

        for (db in dbNames) {
            dbMatch <- object@linksDbId %in% db
            s <- c(s, sprintf("%s: %s", db,
                              paste(IDlinks[dbMatch], IDnames[dbMatch],
                              collapse=", ")))
        }
        cat(strwrap(s, exdent=4), sep="\n")
    })

    ## MATURE
    setClass("mirbaseMatureBimap", contains="AnnDbBimap")
    setClass("mirnaMature",
             representation(matureAccession="character",
                            matureName="character",
                            matureFrom="character",
                            matureTo="character",
                            matureEvidence="character",
                            matureExperiment="character",
                            matureSimilarity="character"))

    setMethod("as.list", "mirbaseMatureBimap", function(x, ...) {
        y <- AnnotationDbi:::flatten(x, fromKeys.only=TRUE)
        makeMBmatureNode <- function(mature_acc, mature_name, mature_from,
                                     mature_to, evidence, experiment,
                                     similarity, ...) {
            new("mirnaMature",
                matureAccession=mature_acc, matureName=mature_name,
                matureFrom=mature_from, matureTo=mature_to,
                matureEvidence=evidence, matureExperiment=experiment,
                matureSimilarity=similarity, ...)
        }
        AnnotationDbi:::.toListOfLists(y, mode=1, makeMBmatureNode)
    })

    setGeneric("matureAccession",
               function(object) standardGeneric("matureAccession"))
    setMethod("matureAccession", "mirnaMature",
              function(object) object@matureAccession)
    setGeneric("matureName", function(object) standardGeneric("matureName"))
    setMethod("matureName", "mirnaMature", function(object) object@matureName)
    ## mature_from and mature_to are stored as char!
    setGeneric("matureFrom", function(object) standardGeneric("matureFrom"))
    setMethod("matureFrom", "mirnaMature",
              function(object) as.integer(object@matureFrom))
    setGeneric("matureTo", function(object) standardGeneric("matureTo"))
    setMethod("matureTo", "mirnaMature",
              function(object) as.integer(object@matureTo))
    setGeneric("matureEvidence",
               function(object) standardGeneric("matureEvidence"))
    setMethod("matureEvidence", "mirnaMature",
              function(object) object@matureEvidence)
    setGeneric("matureExperiment",
               function(object) standardGeneric("matureExperiment"))
    setMethod("matureExperiment", "mirnaMature",
              function(object) object@matureExperiment)
    setGeneric("matureSimilarity",
               function(object) standardGeneric("matureSimilarity"))
    setMethod("matureSimilarity", "mirnaMature",
              function(object) object@matureSimilarity)

    setMethod("show", "mirnaMature", function(object) {
        s <- character(0)
        matureCount <- length(object@matureAccession)
        for (mat in seq(matureCount)) {
            mExp <- object@matureExperiment[mat]
            mSim <- object@matureSimilarity[mat]
            s <- c(s, paste("Accession: ", object@matureAccession[mat], "\n",
                            "  ID: ", object@matureName[mat], "\n",
                            "  Start: ", object@matureFrom[mat], "\n",
                            "  End: ", object@matureTo[mat], "\n",
                            "  Evidence: ", object@matureEvidence[mat], "\n",
                            ifelse(!is.na(mExp),
                                   sprintf("  Experiment: %s\n", mExp), ""),
                            ifelse(mSim != "",
                                   sprintf("  Similarity: %s", mSim), ""),
                            sep=""))
        }
        cat(s, sep="\n")
    })

    ## PMID
    setClass("mirbasePmidBimap", contains="AnnDbBimap")
    setClass("mirnaPmid",
             representation(pmidAuthor="character",
                            pmidTitle="character",
                            pmidJournal="character",
                            pmidMedline="numeric",
                            pmidOrderAdded="numeric"))

    setMethod("as.list", "mirbasePmidBimap", function(x, ...) {
        y <- AnnotationDbi:::flatten(x, fromKeys.only=TRUE)
        makeMBrefNode <- function(medline, author, title, journal,
                                  order_added, ...) {
            new("mirnaPmid",
                pmidAuthor=author, pmidTitle=title,
                pmidJournal=journal, pmidMedline=medline,
                pmidOrderAdded=order_added, ...)
        }
        AnnotationDbi:::.toListOfLists(y, mode=1, makeMBrefNode)
    })

    setGeneric("pmidMedline", function(object) standardGeneric("pmidMedline"))
    setMethod("pmidMedline", "mirnaPmid", function(object) object@pmidMedline)
    setGeneric("pmidAuthor", function(object) standardGeneric("pmidAuthor"))
    setMethod("pmidAuthor", "mirnaPmid", function(object) object@pmidAuthor)
    setGeneric("pmidTitle", function(object) standardGeneric("pmidTitle"))
    setMethod("pmidTitle", "mirnaPmid", function(object) object@pmidTitle)
    setGeneric("pmidJournal", function(object) standardGeneric("pmidJournal"))
    setMethod("pmidJournal", "mirnaPmid", function(object) object@pmidJournal)
    setGeneric("pmidOrderAdded", function(object) {
        standardGeneric("pmidOrderAdded")})
    setMethod("pmidOrderAdded", "mirnaPmid", function(object) {
        object@pmidOrderAdded})

    setMethod("show", "mirnaPmid", function(object) {
        pmidsOrder <- order(object@pmidOrderAdded)
        s <- character(0)
        for (pub in pmidsOrder) {
            sTitle <- object@pmidTitle[pub]
            sPmid <- object@pmidMedline[pub]
            s <- c(s, paste("[", pmidsOrder[pub], "] ",
                            ifelse(sTitle != "", sprintf("%s. ", sTitle), ""),
                            object@pmidAuthor[pub], " ",
                            object@pmidJournal[pub], " ",
                            ifelse(!is.na(sPmid),
                                   sprintf("(PMID:%s)", sPmid), ""), sep=""))
        }
        cat(strwrap(s, exdent=4), sep="\n")
    })


    ## Connect to the SQLite DB
    dbfile <- system.file("extdata", "@DBFILE@", package=pkgname,
                          lib.loc=libname)
    assign("dbfile", dbfile, envir=datacache)
    dbconn <- dbFileConnect(dbfile)
    assign("dbconn", dbconn, envir=datacache)
    ## Create the AnnObj instances
    ann_objs <- createAnnObjs.@DBSCHEMA@("@ANNOBJPREFIX@", "@ANNOBJTARGET@",
                                         dbconn, datacache)
    mergeToNamespaceAndExport(ann_objs, pkgname)
}

.onUnload <- function(libpath) {
    dbFileDisconnect(@ANNOBJPREFIX@_dbconn())
}

MIRBASE_DB_L2Rlink1 <- list(tablename="mirna",
                            Lcolname="mirna_id", Rcolname="_id")

### Mandatory fields: objName, Class and L2Rchain
MIRBASE_DB_AnnDbBimap_seeds <-
    list(
         list(objName="CHR",
              Class="AnnDbBimap",
              L2Rchain=list(
              MIRBASE_DB_L2Rlink1,
              list(tablename="mirna_chromosome_build",
                   Lcolname="_id",
                   Rcolname="xsome"))
              ),
         list(objName="CHRLOC",
              Class="AnnDbBimap",
              L2Rchain=list(
              MIRBASE_DB_L2Rlink1,
              list(tablename="mirna_chromosome_build",
                   Lcolname="_id",
                   tagname=c(chromosome="{xsome}"),
                   Rcolname="contig_start"))
              ),
         list(objName="CHRLOCEND",
              Class="AnnDbBimap",
              L2Rchain=list(
              MIRBASE_DB_L2Rlink1,
              list(tablename="mirna_chromosome_build",
                   Lcolname="_id",
                   tagname=c(chromosome="{xsome}"),
                   Rcolname="contig_end"))
              ),
         list(objName="CLUSTER",
              Class="AnnDbBimap",
              L2Rchain=list(
              list(tablename="mirna_cluster",
                   Lcolname="mirna_id",
                   Rcolname="member"))
              ),
         list(objName="COMMENT",
              Class="AnnDbBimap",
              L2Rchain=list(
              list(tablename="mirna",
                   Lcolname="mirna_id",
                   Rcolname="comment"))
              ),
         list(objName="CONTEXT",
              Class="mirbaseContextBimap",
              L2Rchain=list(
              MIRBASE_DB_L2Rlink1,
              list(tablename="mirna_context",
                   Lcolname="_id",
                   Rcolname="transcript_id",
                   Rattribnames=c(overlap_sense="{overlap_sense}",
                   overlap_type="{overlap_type}",
                   number="{number}",
                   transcript_source="{transcript_source}",
                   transcript_name="{transcript_name}")))
              ),
         list(objName="DESCRIPTION",
              Class="AnnDbBimap",
              L2Rchain=list(
              list(tablename="mirna",
                   Lcolname="mirna_id",
                   Rcolname="description"))
              ),
         list(objName="FAMILY",
              Class="AnnDbBimap",
              L2Rchain=list(
              MIRBASE_DB_L2Rlink1,
              list(tablename="mirna_2_prefam",
                   Lcolname="_id",
                   Rcolname="auto_prefam"),
              list(tablename="mirna_prefam",
                   Lcolname="auto_prefam",
                   tagname=c(id="{prefam_id}"),
                   Rcolname="prefam_acc"))
              ),
         list(objName="HAIRPIN",
              Class="AnnDbBimap",
              L2Rchain=list(
              list(tablename="mirna_hairpin",
                   Lcolname="mirna_id",
                   Rcolname="hairpin"))
              ),
         list(objName="ID2ACC",
              Class="AnnDbBimap",
              L2Rchain=list(
              list(tablename="mirna",
                   Lcolname="mirna_id",
                   Rcolname="mirna_acc"))
              ),
         list(objName="ID2SPECIES",
              Class="AnnDbBimap",
              L2Rchain=list(
              list(tablename="mirna",
                   Lcolname="mirna_id",
                   Rcolname="organism"))
              ),
         list(objName="LINKS",
              Class="mirbaseLinksBimap",
              L2Rchain=list(
              MIRBASE_DB_L2Rlink1,
              list(tablename="mirna_database_links",
                   Lcolname="_id",
                   tagname=c(db_id="{db_id}"),
                   Rcolname="db_link",
                   Rattribnames=c(db_secondary="{db_secondary}")))
              ),
         list(objName="MATURE",
              Class="mirbaseMatureBimap",
              L2Rchain=list(
              MIRBASE_DB_L2Rlink1,
              list(tablename="mirna_pre_mature",
                   Lcolname="_id",
                   Rcolname="auto_mature"),
              list(tablename="mirna_mature",
                   Lcolname="auto_mature",
                   tagname=c(mature_name="{mature_name}"),
                   Rcolname="mature_acc",
                   Rattribnames=c(mature_from="{mature_from}",
                   mature_to="{mature_to}",
                   evidence="{evidence}",
                   experiment="{experiment}",
                   similarity="{similarity}")))
              ),
         list(objName="MFE",
              Class="AnnDbBimap",
              L2Rchain=list(
              list(tablename="mirna_hairpin",
                   Lcolname="mirna_id",
                   Rcolname="mfe"))
              ),
         list(objName="PMID",
              Class="mirbasePmidBimap",
              L2Rchain=list(
              MIRBASE_DB_L2Rlink1,
              list(tablename="mirna_literature_references",
                   Lcolname="_id",
                   Rcolname="order_added",
                   Rattribnames=c(medline="literature_references.medline",
                   author="literature_references.author",
                   title="literature_references.title",
                   journal="literature_references.journal"),
                   Rattrib_join=paste("LEFT JOIN literature_references ON",
                   "{auto_lit}=literature_references.auto_lit", sep=" ")))
              ),
         list(objName="SEQUENCE",
              Class="AnnDbBimap",
              L2Rchain=list(
              MIRBASE_DB_L2Rlink1,
              list(tablename="mirna",
                   Lcolname="_id",
                   Rcolname="sequence"))
              ),
         list(objName="SPECIES",
              Class="AnnDbBimap",
              L2Rchain=list(
              list(tablename="mirna_species",
                   Lcolname="organism",
                   Rcolname="name",
                   Rattribnames=c(division="{division}",
                   taxonomy="{taxonomy}",
                   genome_assembly="{genome_assembly}",
                   ensembl_db="{ensembl_db}"))
              )
         )
)

createAnnObjs.MIRBASE_DB <- function(prefix, objTarget, dbconn, datacache) {
    AnnotationDbi:::checkDBSCHEMA(dbconn, "MIRBASE_DB")

    ## AnnDbBimap objects
    seed0 <- list(objTarget=objTarget, datacache=datacache)
    ann_objs <- AnnotationDbi:::createAnnDbBimaps(MIRBASE_DB_AnnDbBimap_seeds,
                                                  seed0)

    ## Reverse maps
    ann_objs$ACC2ID <- AnnotationDbi:::revmap(ann_objs$ID2ACC, objName="ACC2ID")
    ann_objs$SPECIES2ID <- AnnotationDbi:::revmap(ann_objs$ID2SPECIES,
                                                  objName="SPECIESID")

    ## 1 special map that is not an AnnDbBimap object
    ##(just a named integer vector)
    ann_objs$MAPCOUNTS <- AnnotationDbi:::createMAPCOUNTS(dbconn, prefix)

    AnnotationDbi:::prefixAnnObjNames(ann_objs, prefix)
}



