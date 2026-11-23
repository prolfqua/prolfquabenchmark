#!/usr/bin/env Rscript
# Writes the input of time_prolfqua_steps.R: the IonStar MaxQuant peptides joined with the sample annotation, read
# with the installed prolfquapp. Prepared once, so the timed prolfqua versions all start from the same table.
#   Rscript workflow/scripts/timing_input.R <out.rds>

out <- commandArgs(trailingOnly = TRUE)[[1]]
datadir <- file.path(find.package("prolfquadata"), "quantdata")
peptides <- prolfquapp::tidyMQ_Peptides(file.path(datadir, "MAXQuant_IonStar2018_PXD003881.zip"))
annotation <- readxl::read_xlsx(file.path(datadir, "annotation_Ionstar2018_PXD003881.xlsx"))
saveRDS(dplyr::inner_join(peptides, annotation, by = "raw.file"), out)
