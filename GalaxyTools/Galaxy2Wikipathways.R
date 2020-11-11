if(!"rWikiPathways" %in% installed.packages()){
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("rWikiPathways", update = FALSE)
}

load.libs <- c(
  "RColorBrewer",
  "rWikiPathways",
  "RCy3",
  "RCurl")

options(install.packages.check.source = "no")
options(install.packages.compile.from.source = "never")
if (!require("pacman")) install.packages("pacman"); library(pacman)
p_load(load.libs, update = TRUE, character.only = TRUE)
status <- sapply(load.libs,require,character.only = TRUE)
if(!all(status)){
  cat("ERROR: One or more libraries failed to install correctly. Check the following list for FALSE cases and try again...\n\n")
  status
}

# inputs
wp_id <- 'WP528'
headers <- TRUE
data_url <- 'http://localhost:3000/galaxy_test_Copy.txt'

#fetch and modify data
dat <- read.table(text=getURL(data_url), header = headers, fill=TRUE)
if(headers)
{
  colnames(dat)[1] <- "geneid"
  colnames(dat)[2] <- "fc"
} else
{
  colnames(dat) <- c("geneid", "fc", "p-value", "adj p-value")
}

dat <- data.frame(lapply(dat, function(v) {
  if (is.character(v)) return(toupper(v))
  else return(v)
}))
min.fc = min(dat["fc"],na.rm=TRUE)
max.fc = max(dat["fc"],na.rm=TRUE)
abs.fc = max(abs(min.fc),abs(max.fc))
data.values = c(-abs.fc,0,abs.fc)

#open in browser
rgb2hex <- function(r,g,b) rgb(r, g, b, maxColorValue = 255)
fill_hex <- function(x) {
  out <- vector(mode="character", length=length(x))
  for (i in seq_along(x)) 
  {
    if(is.na(x[[i]]))
    {
      out[[i]] <- "#ffffff"
    }
    else
    {
      rel_value = 70+(185*( 1 - (abs(x[[i]])/abs.fc)))
      if(x[[i]] < 0)
      {
        if(abs(min.fc) == 0)
        {
          out[[i]] <- "#ffffff"
        }
        else{
          out[[i]] <- rgb2hex(rel_value, rel_value, 255)
        }
      }
      else
      {
        if(abs(min.fc) == 0)
        {
          out[[i]] <- "#ffffff"
        }
        else{
          out[[i]] <- rgb2hex(255, rel_value, rel_value)
        }
      }
    }
  

  }
  return(out)
}
dat$hex <- fill_hex(dat[["fc"]])
url_params <- paste(with(dat, paste(substr(hex,2,nchar(hex)), geneid, sep="=HGNC_")), collapse="&")
browseURL(paste(c("https://pathway-viewer.toolforge.org/?id=", wp_id, "&", url_params), collapse = ""), browser = getOption("browser"))

#open in CytoScape
tryCatch(
{
  cytoscapePing()
  installApp('WikiPathways') 
  RCy3::commandsRun(paste(c('wikipathways import-as-pathway id=',wp_id), collapse = ""))
  toggleGraphicsDetails()
  node.colors <- c(rev(brewer.pal(length(data.values), "RdBu")))
  loadTableData(dat, table = "node", data.key.column = "geneid", table.key.column = "name")
  setNodeColorMapping("fc", data.values, node.colors, default.color = "#FFFFFF", style.name = "WikiPathways")
},
error=function(cond) {
  # catching if CytoScape is not launched
})
