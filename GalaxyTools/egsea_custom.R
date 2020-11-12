library(EGSEA)

genes <- grep("^#", readLines("input.txt"), invert = T, value = T)

library(jsonlite)
dmaps_out <- fromJSON(paste(readLines("dmaps_out.json"), collapse = ""))
dmaps_list <- sapply(dmaps_out[[2]], function(x) x[["hgncs"]])
names(dmaps_list) <- sapply(dmaps_out[[2]], function(x) x[["name"]])
### Build a custom index for the disease map gene sets
dmap.gs.annot <- buildCustomIdx(geneIDs = genes, gsets = dmaps_list, species="human")
### Run egsea.ora for set-based enrichment (overrepresentation analysis), with report = F to grab the results
egsea.result <- egsea.ora(geneIDs = genes, universe = unique(unlist(dmaps_list)), gs.annots = dmap.gs.annot, report = F)
### Parse the results
### 1) Grab the Experimental Contrast part of the result
ec <- egsea.result@results$custom$test.results$ExperimentalContrast

names



wps_out <- fromJSON(paste(readLines("wikipathways.json"), collapse = ""))
wps_list <- wps_out[[2]]
names(wps_list) <- unlist(wps_out[[1]])
wp.gs.annot = buildCustomIdx(geneIDs = genes, gsets = wps_list, species="human")
egsea.ora(geneIDs = genes, universe = unique(unlist(wps_list)), gs.annots = wp.gs.annot)
