### Map and pathway enrichment analysis.
### Requires commandline input or file "input.txt" to be present in the script directory.

required.packages <- c("jsonlite", "httr")
new.packages <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos='http://cran.us.r-project.org')

library(jsonlite)
library(httr)

### A convenience function to fetch data from MINERVA
ask_GET <- function(furl, fask) {
  resp <- httr::GET(url = paste0(furl, fask),
                    httr::add_headers('Content-Type' = "application/x-www-form-urlencoded"))
  if(httr::status_code(resp) == 200) {
    return(httr::content(resp, as = "text"))
  }
  return(NULL)
}

### In case we want to pass parameters as an argument
args = commandArgs(trailingOnly=TRUE)

### Config listing the resources, which we want to fetch
config <- "config.txt"
if(length(args) > 0) { config <- args[1] }
if(!file.exists(config)) { message("No correct config given!")}
config <- read.table(config, sep = ",", stringsAsFactors = F, header = T)

### for a given model id (for all maps/submaps in the project)
### Fetch HGNC symbols for proteins, RNAs and genes.
get_hgncs <- function(fmodelid, all_models) {
  ### MINERVA call fetching elements with relevant parameters
  fmodel <- fromJSON(ask_GET(paste0(mnv_base,"models/",fmodelid,"/"), "bioEntities/elements/?columns=id,name,type,references,bounds"), flatten = F)
  
  ### Model name used in naming gene sets
  model_name <- all_models[all_models$idObject == fmodelid,"name"]
  
  ### If no annotations, break
  if(all(is.null(unlist(fmodel$references)))) { return(NULL) }
  
  ### List all HGNCs in maps and submaps
  ret <- list(
    list(name = model_name, 
         hgncs = unlist(sapply(fmodel$references[fmodel$type %in% c("Protein", "RNA", "Gene")], function(x) x[x$type == "HGNC_SYMBOL", "resource"]))))
  
  ### List all HGNCs in pathways (smaller areas within a single diagram)
  if(length(unique(fmodel[fmodel$type == "Pathway","name"])) > 0) {
    ### Calculating the inclusion of the elements per bioentity set
    for(pname in unique(me[fmodel$type == "Pathway","name"])) {
      res <- apply(me[fmodel$name == pname,"bounds"], 1, 
                   function(x) fmodel$bounds$x >= x["x"] & fmodel$bounds$x <= (x["x"] +  x["width"]) & 
                     fmodel$bounds$y >= x["y"] & fmodel$bounds$y <= (x["y"] +  x["height"]))
      
      hgncs <- unlist(sapply(fmodel$references[res & fmodel$type %in% c("Protein", "RNA", "Gene")], 
                             function(x) x[x$type == "HGNC_SYMBOL", "resource"]))
      ret <- c(ret, list(name = paste0(model_name,":",pname), hgncs = hgncs))
    }
  }
  return(ret)
}

message("Fetching contents of disease maps in config")

all_maps <- list()

for(map in config[config$type == "map","resource"]) {
  message(paste0("Querying ", map))
  cfg <- ask_GET(map, "configuration/")
  ### For erroneous response, skip to next map
  if(is.null(cfg)) { next }
  cfg <- fromJSON(cfg)
  project_id <- cfg$options[cfg$options$type == "DEFAULT_MAP","value"]
  
  mnv_base <- paste0(map,"projects/",project_id,"/")
  
  ### Ask for models
  message(paste0("Asking for models: ", mnv_base, "models/"))
  models <- ask_GET(mnv_base, "models/")
  ### For erroneous response, skip to next map
  if(is.null(models)) { next }
  models <- fromJSON(models, flatten = F)
  
  message(paste0("Fetching HGNCs for individual models..."))
  
  model_hgncs <- lapply(models$idObject, get_hgncs, models)
  
  #Remove null entries
  model_hgncs <- model_hgncs[!sapply(model_hgncs, is.null)]
 
  all_maps <- c(all_maps, list(map_name = map, map_elements = model_hgncs))
   
  message("Done.")
}

message("Writing disease map(s) to JSON...")
cat(toJSON(all_maps), file = "dmaps_out.json")

### Massive SPARQL query (https://www.wikipathways.org/index.php/Help:WikiPathways_Sparql_queries#List_of_WikiPathways_for_HGNC_symbols)
### with the Homo sapiens condition added

sparql_url <- paste0("http://sparql.wikipathways.org/sparql?&query=select+distinct+%3FpathwayRes+str%28%3Fwpid%29+as+%3Fpathway+str%28%3Ftitle%29+as+%3FpathwayTitle+fn%3Asubstring%28%3FhgncId%2C36%29+as+%3FHGNC+where+%7B%0D%0A",
                     "++%3Fgene+a+wp%3AGeneProduct+%3B%0D%0A",
                     "++++dcterms%3Aidentifier+%3Fid+%3B%0D%0A",
                     "++++dcterms%3AisPartOf+%3FpathwayRes+%3B%0D%0A",
                     "++++wp%3AbdbHgncSymbol+%3FhgncId+.%0D%0A",
                     "++%3FpathwayRes+a+wp%3APathway+%3B%0D%0A",
                     "++++wp%3AorganismName+%22Homo%20sapiens%22%5E%5Exsd%3Astring+%3B%0D%0A", #Homo sapiens condition here
                     "++++dcterms%3Aidentifier+%3Fwpid+%3B%0D%0A",
                     "++++dc%3Atitle+%3Ftitle+.%0D%0A%7D",
                     "&format=text%2Fcsv&timeout=0&debug=on")

### Retrieve the table whole opening/closing the connection gracefully
con <- url(sparql_url)
res <- try(read.table(con, sep = ",", header = T, stringsAsFactors = F), 
               silent = T)
if(class(res) == "try-error") {
  message(paste0("Cannot read from ", sparql_url))
  message(res)
  close(con)
}

### Compile a similar JSON with the WikiPathways pathways

### Melt the table by the unique pathway name
wp_list <- lapply(unique(res$pathwayTitle), function(x) list(name = x, hgncs = unique(res[res$pathwayTitle == x, "HGNC"])))

message("Writing WikiPathways to JSON...")
cat(toJSON(wp_list), file = "wikipathways.json")

message("Done.")


### Spare code ###

### SPARQL query using R library; not optimal, the library is from 2013

# library(SPARQL)
# 
# endpoint <- "http://sparql.wikipathways.org/sparql"
# 
# q <- "select distinct ?pathwayRes (str(?wpid) as ?pathway) (str(?title) as ?pathwayTitle) (fn:substring(?hgncId,36) as ?HGNC) where {
#           ?gene a wp:GeneProduct ;
#             dcterms:identifier ?id ;
#             dcterms:isPartOf ?pathwayRes ;
#             wp:bdbHgncSymbol ?hgncId .
#           ?pathwayRes a wp:Pathway ;
#             dcterms:identifier ?wpid ;
#             dc:title ?title .
#           ?pathway wp:organismName 'Homo sapiens^^xsd:string .
#       }"
# 
# prefix <- c("wp","http://vocabularies.wikipathways.org/wp#",
#             "rdfs","http://www.w3.org/2000/01/rdf-schema#",
#             "dcterms","http://purl.org/dc/terms/",
#             "xsd","http://www.w3.org/2001/XMLSchema#")
# 
# res <- SPARQL(endpoint,q,prefix)$results
# 
# 
