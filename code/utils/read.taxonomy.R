read.taxonomy <- function(file) {
  dt <- data.table::fread(file, header = FALSE, sep = "\t", col.names = c("otu", "taxon"))
  dt[, c("kingdom", "kingdom.pcnt", "phylum", "phylum.pct", "class", "class.pct",
         "order", "order.pct", "family", "family.pct", "genus", "genus.pct") :=
       tstrsplit(taxon, "\\(|;|\\)\\;")]
  dt[, taxon := NULL]
  dt
}
