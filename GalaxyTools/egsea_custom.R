library(BiocManager)

BiocManager::install("EGSEA")

library(EGSEA)

genes <- grep("^#", readLines("input.txt"), invert = T, value = T)

library(jsonlite)
dmaps_out <- fromJSON(paste(readLines("dmaps_out.json"), collapse = ""))
dmaps_list <- sapply(dmaps_out[[2]], function(x) x[["hgncs"]])
names(dmaps_list) <- sapply(dmaps_out[[2]], function(x) x[["name"]])
dmap.gs.annot = buildCustomIdx(geneIDs = genes, gsets = dmaps_list, species="human")
egsea.ora(geneIDs = genes, universe = unique(unlist(dmaps_list)), gs.annots = dmap.gs.annot)

wps_out <- fromJSON(paste(readLines("wikipathways.json"), collapse = ""))
wps_list <- wps_out[[2]]
names(wps_list) <- unlist(wps_out[[1]])
wp.gs.annot = buildCustomIdx(geneIDs = genes, gsets = wps_list, species="human")
egsea.ora(geneIDs = genes, universe = unique(unlist(wps_list)), gs.annots = wp.gs.annot)
