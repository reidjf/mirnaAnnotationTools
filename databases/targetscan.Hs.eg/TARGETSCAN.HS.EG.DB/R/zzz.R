datacache <- new.env(hash=TRUE, parent=emptyenv())

## append targetscan version and date
@ANNOBJPREFIX@ <- function() {
    showQCData("@ANNOBJPREFIX@", datacache)
    cat(paste("targtscan version: ",
              dbmeta(datacache, 'TARGETSCAN.HS.EGSOURCEVERSION'), "\n", sep=""))
    cat(paste("targetscan date: ",
              dbmeta(datacache, 'TARGETSCAN.HS.EGSOURCEDATE'), "\n", sep=""))
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

    setClass("miRNATargetAnnDbBimap", contains="AnnDbBimap")

    setClass("targetscanTarget",
             representation(miRFamily="character",
                            UTRstart="numeric",
                            UTRend="numeric",
                            MSAstart="numeric",
                            MSAend="numeric",
                            Seedmatch="character",
                            PCT="character"))

    setMethod("as.list", "miRNATargetAnnDbBimap", function(x, ...) {
        y <- AnnotationDbi:::flatten(x, fromKeys.only=TRUE)
        makemiRNATargetNode <- function(name, UTR_start, UTR_end,
                                        MSA_start, MSA_end, Seed_match,
                                        PCT, ...) {
            new("targetscanTarget",
                miRFamily=name, UTRstart=UTR_start, UTRend=UTR_end,
                MSAstart=MSA_start, MSAend=MSA_end, Seedmatch=Seed_match,
                PCT=PCT)
        }
        AnnotationDbi:::.toListOfLists(y, mode=1, makemiRNATargetNode)
    })

    setGeneric("miRFamily", function(object) standardGeneric("miRFamily"))
    setMethod("miRFamily", "targetscanTarget",
              function(object) object@miRFamily)
    setGeneric("UTRstart", function(object) standardGeneric("UTRstart"))
    setMethod("UTRstart", "targetscanTarget", function(object) object@UTRstart)
    setGeneric("UTRend", function(object) standardGeneric("UTRend"))
    setMethod("UTRend", "targetscanTarget", function(object) object@UTRend)
    setGeneric("MSAstart", function(object) standardGeneric("MSAstart"))
    setMethod("MSAstart", "targetscanTarget", function(object) object@MSAstart)
    setGeneric("MSAend", function(object) standardGeneric("MSAend"))
    setMethod("MSAend", "targetscanTarget", function(object) object@MSAend)
    setGeneric("Seedmatch", function(object) standardGeneric("Seedmatch"))
    setMethod("Seedmatch", "targetscanTarget",
              function(object) object@Seedmatch)
    setGeneric("PCT", function(object) standardGeneric("PCT"))
    setMethod("PCT", "targetscanTarget", function(object) object@PCT)

    setMethod("show", "targetscanTarget", function(object) {
        tT <- length(object@miRFamily)
        s <- character(0)
        for (tG in 1:tT) {
            for (slotname in slotNames(object)) {
                x <- slot(object, slotname)[tG]
                if (length(x) == 0)  next
                s <- c(s, paste(slotname, ": ", x, sep=""))
            }
            s <- c(s, "\n")
        }
        cat(strwrap(s, exdent=4), sep="\n")
    })

    setClass("miRNAAnnDbBimap", contains="AnnDbBimap")

    setClass("targetscanMiRBase",
             representation(MiRBaseID="character",
                            MiRBaseAccession="character",
                            Seedm8="character",
                            Species="character",
                            Maturesequence="character",
                            Familyconservation="character"))

    setMethod("as.list", "miRNAAnnDbBimap", function(x, ...) {
        y <- AnnotationDbi:::flatten(x, fromKeys.only=TRUE)
        makemiRNANode <- function(mirbase_id, MiRBase_Accession,
                                  Seed_m8, Species_ID, Mature_sequence,
                                  Family_conservation, ...) {
            new("targetscanMiRBase",
                MiRBaseID=mirbase_id, MiRBaseAccession=MiRBase_Accession,
                Seedm8=Seed_m8, Species=Species_ID,
                Maturesequence=Mature_sequence,
                Familyconservation=Family_conservation)
        }
        AnnotationDbi:::.toListOfLists(y, mode=1, makemiRNANode)
    })

    setGeneric("MiRBaseID", function(object) standardGeneric("MiRBaseID"))
    setMethod("MiRBaseID", "targetscanMiRBase",
              function(object) object@MiRBaseID)
    setGeneric("MiRBaseAccession", function(object)
               standardGeneric("MiRBaseAccession"))
    setMethod("MiRBaseAccession", "targetscanMiRBase",
              function(object) object@MiRBaseAccession)
    setGeneric("Seedm8", function(object) standardGeneric("Seedm8"))
    setMethod("Seedm8", "targetscanMiRBase",
              function(object) object@Seedm8)
    setGeneric("Species", function(object) standardGeneric("Species"))
    setMethod("Species", "targetscanMiRBase",
              function(object) object@Species)
    setGeneric("Maturesequence",
               function(object) standardGeneric("Maturesequence"))
    setMethod("Maturesequence", "targetscanMiRBase",
              function(object) object@Maturesequence)
    setGeneric("Familyconservation",
               function(object) standardGeneric("Familyconservation"))
    setMethod("Familyconservation", "targetscanMiRBase",
              function(object) object@Familyconservation)

    setMethod("show", "targetscanMiRBase", function(object) {
        for (slotname in slotNames(object)) {
            name2 <- sub("\\.", " ", slotname)
            s <- paste(sep="", name2, ": ", slot(object, slotname))
            cat(strwrap(s, exdent=4), sep="\n")
        }
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

### Mandatory fields: objName, Class and L2Rchain
TARGETSCAN_DB_AnnDbBimap_seeds <- list(
    list(objName="MIRBASE2FAMILY",
         Class="AnnDbBimap",
         L2Rchain=list(
           list(tablename="mirbase",
                Lcolname="mirbase_id",
                Rcolname="family"
                ),
           list(
                tablename="mirna_family",
                Lcolname="_id",
                Rcolname="name"
                )
           )
         ),
     list(objName="MIRNA",
          Class="miRNAAnnDbBimap",
          L2Rchain=list(
            list(tablename="mirbase",
                 Lcolname="mirbase_id",
                 Rcolname="mirbase_id",
                 Rattribnames=c(
                   MiRBase_Accession="{mirbase_accession}",
                   Seed_m8="{seed_m8}",
                   Species_ID="species.name",
                   Mature_sequence="{mature_sequence}",
                   Family_conservation="{family_conservation}"),
                 Rattrib_join="LEFT JOIN species ON {species}=species.id"
                 )
            )
          ),
      list(objName="TARGETS",
           Class="AnnDbBimap",
           L2Rchain=list(
             list(tablename="genes",
                  Lcolname="gene_id",
                  Rcolname="_id"
                  ),
             list(tablename="targets",
                  Lcolname="target",
                  Rcolname="family"
                  ),
             list(tablename="mirna_family",
                  Lcolname="_id",
                  Rcolname="name"
                  )
             )
           ),
      list(objName="TARGETSFULL",
           Class="miRNATargetAnnDbBimap",
           L2Rchain=list(
             list(tablename="genes",
                  Lcolname="gene_id",
                  Rcolname="_id"
                  ),
             list(tablename="targets",
                  Lcolname="target",
                  Rattribnames=c(
                    UTR_start="{utr_start}",
                    UTR_end="{utr_end}",
                    MSA_start="{msa_start}",
                    MSA_end="{msa_end}",
                    Seed_match="seed_match.name",
                    PCT="{pct}"),
                  Rattrib_join="LEFT JOIN seed_match ON {seed_match}=seed_match._id LEFT JOIN mirna_family AS _R ON {family}=_R._id",
                  Rcolname="name"
##                   ),
##              list(tablename="mirna_family",
##                   Lcolname="_id",
##                   Rcolname="name"
                  )
             )
           )
)

createAnnObjs.TARGETSCAN.HS.EG_DB <- function(prefix, objTarget,
                                               dbconn, datacache) {
    ## checkDBSCHEMA(dbconn, "TARGETSCAN_DB")

    ## AnnDbBimap objects
    seed0 <- list(
                  objTarget=objTarget,
                  datacache=datacache
                  )
    ann_objs <- AnnotationDbi:::createAnnDbBimaps(TARGETSCAN_DB_AnnDbBimap_seeds, seed0)

    ## Reverse maps
    revmap2 <- function(from, to)
    {
        map <- revmap(ann_objs[[from]], objName=to)
        L2Rchain <- map@L2Rchain
        tmp <- L2Rchain[[1]]@filter
        L2Rchain[[1]]@filter <-
          L2Rchain[[length(L2Rchain)]]@filter
        L2Rchain[[length(L2Rchain)]]@filter <- tmp
        map@L2Rchain <- L2Rchain
        map
    }
    ann_objs$FAMILY2MIRBASE <- revmap2("MIRBASE2FAMILY", "FAMILY2MIRBASE")

    ## 1 special map that is not an AnnDbBimap object (just a named integer vector)
    ann_objs$MAPCOUNTS <- AnnotationDbi:::createMAPCOUNTS(dbconn, prefix)

    AnnotationDbi:::prefixAnnObjNames(ann_objs, prefix)
}

