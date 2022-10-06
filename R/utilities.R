tolong <- function(assay, rowData, colData ){
    annot <- as_tibble(data.frame(colData))
    hierarchy <- prolfqua::matrix_to_tibble( data.frame( rowData ) )
    ll <- prolfqua::matrix_to_tibble(assay)
    ll <- ll |> pivot_longer(cols = -all_of("row.names"), names_to = "sampleName", values_to = "intensity")
    hl <- inner_join(hierarchy, ll)
    ahl <- inner_join(annot, hl)
}
