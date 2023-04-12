
#' get all known M.oryzae gene ids - the universe
#' @param genefile csv file with GeneID column of all Mo ids
all_mo_genes <- function(genefile=NA){
  #here::here("inst", "extdata", "all_mo_ids.csv")
  if (! is.na(genefile)){
    readr::read_csv(genefile)$GeneID
  }
  else {
    mo_gene_ids
  }
}

#' get mapping between Mo Gene IDs and GO Terms
#' @param termsfile TSV file from biomart export of GO terms
mo_mapping <- function(termsfile=NA) {
  #here::here("inst", "extdata", "MG8_mart_export.txt")
  terms <- terms_data
  if (! is.na(termsfile)){
    terms <- readr::read_tsv(termsfile)
  }
  term2gene <- data.frame(
    term = terms$`GO term accession`,
    gene = terms$`Gene stable ID`
  )

  term2name <- data.frame(
    term = terms$`GO term accession`,
    name = terms$`GO term name`
  )

  all_genes <- all_mo_genes()

  return(list(
    term2gene = term2gene,
    term2name = term2name,
    all_genes = all_genes
  ))
}

#' run clusterProfiler::enricher on vector of gene ids,
#' @param genes character vector of gene IDs of interest
#' @param termsfile TSV file from biomart export of GO terms
#' @return enricher object
#' @export
do_enrich <- function(genes, termsfile=NA, ...) {

  here::here("inst", "extdata", "MG8_mart_export.txt")
  info <- mo_mapping()
  if (! is.na(termsfile)){
    info <- mo_mapping(termsfile=termsfile)
  }
  clusterProfiler::enricher( genes,
                             universe=info$all_genes,
                             TERM2GENE=info$term2gene,
                             TERM2NAME=info$term2name,
                             ...
                             )

}

#' conver enricher object to DAVID format
#' @param enrich enricher object from do_enrich
#' @param termsfile TSV file from biomart export of GO terms
#' @return DAVID format dataframe
#' @export
enricher_to_david <- function(enrich,termsfile=NA){
  comma_sep_genes = gsub("/", ", ", enrich@result$geneID)
  #here::here("inst","extdata", "MG8_mart_export.txt")
  terms <- terms_data
  if (!is.na(termsfile)){
    terms <- readr::read_tsv(termsfile)

  terms <- dplyr::mutate(terms,
    short_category = dplyr::if_else(`GO domain` == 'biological_process', "BP",
                                    dplyr::if_else(`GO domain` == "molecular_function", "MF", "CC"))
  )

  id_to_category <- terms$short_category
  names(id_to_category) <- terms$`GO term accession`
  fixed_category <- id_to_category[enrich@result$ID]

  data.frame(
    category = fixed_category,
    ID = enrich@result$ID,
    term = enrich@result$Description,
    genes = comma_sep_genes,
    adj_pval = enrich@result$p.adjust
  )


}
