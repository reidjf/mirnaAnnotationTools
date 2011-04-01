library("AnnotationDbi")

## this is the devel-version of makeAnnDbPkg (1.9.4)
## in '2.5' release we have
## makeAnnDbPkg(seed, dbfile, dest_dir = tempdir())
## Error in createPackage(x@Package, destinationDir = dest_dir, originDir = template_path,  :
##  'originDir' must be a directory ()
makeAnnDbPkg <- function (x, dbfile, dest_dir = ".", no.man = FALSE, ...) {
    x <- AnnotationDbi:::initWithDbMetada(x, dbfile)
    x <- AnnotationDbi:::initComputedSlots(x)
    dbfile_basename <- basename(dbfile)
    if (dbfile_basename != paste(x@AnnObjPrefix, ".sqlite", sep = ""))
        stop("'", dbfile, "': File name doesn't match 'x@AnnObjPrefix' (",
             x@AnnObjPrefix, ")")
##jfr
    ##if (!grepl("^/", x@PkgTemplate)[1]) {
    ##    template_path <- system.file("AnnDbPkg-templates",
    ##                                 x@PkgTemplate,
    ##                                 package="AnnotationDbi")
    ##} else {
        template_path <- x@PkgTemplate
    ##}

    ann_dbi_version <- installed.packages()["AnnotationDbi", "Version"]
    symvals <- list(
                   DBSCHEMA = x@DBschema,
                   PKGTITLE = x@Title,
                   ANNOBJPREFIX = x@AnnObjPrefix,
                   ANNOBJTARGET = x@AnnObjTarget,
                   ORGANISM = x@organism,
                   SPECIES = x@species,
                   MANUF = x@manufacturer,
                   CHIPNAME = x@chipName,
                   MANUFURL = x@manufacturerUrl,
                   AUTHOR = x@Author,
                   MAINTAINER = x@Maintainer,
                   PKGVERSION = x@Version,
                   LIC = x@License,
                   BIOCVIEWS = x@biocViews,
                   DBFILE = dbfile_basename,
                   ANNDBIVERSION = ann_dbi_version
                    )
    man_dir <- file.path(template_path, "man")
    if (file.exists(man_dir)) {
        if (!no.man) {
            doc_template_names <- list.files(man_dir, "\\.Rd$")
            #is_static <- doc_template_names %in% c("_dbconn.Rd", "_dbfile.Rd")
            #doc_template_names <- doc_template_names[!is_static]
            map_names <- sub("\\.Rd$", "", doc_template_names)
            if (length(map_names) != 0)
                symvals <- c(symvals,
                             AnnotationDbi:::getSymbolValuesForManPages(map_names, dbfile))
        }
        else {
            unlink(man_dir, recursive = TRUE)
        }
    }
    if (any(duplicated(names(symvals)))) {
        str(symvals)
        stop("'symvals' contains duplicated symbols (see above)")
    }
    Biobase::createPackage(x@Package,
                           destinationDir = dest_dir,
                           originDir = template_path,
                           symbolValues = symvals)
    ## rename Rd files
    if (file.exists(man_dir) && !no.man && length(doc_template_names) != 0) {
        doc_path <- file.path(dest_dir, x@Package, "man")
        from_doc_names <- paste(doc_path, doc_template_names, sep = .Platform$file.sep)
        to_doc_names <- paste(x@AnnObjPrefix, doc_template_names, sep = "")
        to_doc_names <- paste(doc_path, to_doc_names, sep = .Platform$file.sep)
        mapply(file.rename, from_doc_names, to_doc_names)
    }

    dest_db_dir <- file.path(dest_dir, x@Package, "inst", "extdata")
    if (!file.exists(dest_db_dir) &&
        !dir.create(dest_db_dir, recursive = TRUE))
        stop("unable to create dest db dir ", dest_db_dir)
    dest_dbfile <- file.path(dest_db_dir, dbfile_basename)
    if (!file.copy(dbfile, dest_dbfile))
        stop("cannot copy file '", dbfile, "' to '", dest_dbfile, "'")
    if (.Platform$OS.type != "windows") {
        command <- paste("chmod 444", dest_dbfile)
        if (system(command) != 0)
            warning(command, " failed")
    }
    return(invisible(TRUE))
}

