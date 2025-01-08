tolong <- function(assay, rowData, colData ){
    annot <- as_tibble(data.frame(colData))
    hierarchy <- prolfqua::matrix_to_tibble( data.frame( rowData ) )
    ll <- prolfqua::matrix_to_tibble(assay)
    ll <- ll |> pivot_longer(cols = -all_of("row.names"), names_to = "sampleName", values_to = "intensity")
    hl <- inner_join(hierarchy, ll)
    ahl <- inner_join(annot, hl)
}


#' preprocess cptac dataset for benchmarking
#' @export
cptac_bench_preprocess <- function(data, idcol = "protein_Id") {
    tmp <- data |>
        dplyr::ungroup() |>
        dplyr::mutate(species  = dplyr::case_when(
            grepl("YEAST", !!sym(idcol)) ~ "YEAST",
            grepl("UPS", !!sym(idcol)) ~ "UPS",
            TRUE ~ "OTHER"
        ))
    res <- tmp |> dplyr::filter(!.data$species == "OTHER")
    res <- res |> dplyr::mutate(TP = (.data$species == "UPS"))
    return(list(data = res , table = table(tmp$species)))
}
