library(curl)

download.file("https://gitlab.bfabric.org/wolski/prolfquadata/-/raw/master/inst/quantdata/MAXQuant_IonStar2018_PXD003881.zip?inline=false", "MAXQuant_IonStar2018_PXD003881.zip")
tmp <- unzip("MAXQuant_IonStar2018_PXD003881.zip", list = RUE)
pg <- read.table( unz("MAXQuant_IonStar2018_PXD003881.zip", "proteinGroups.txt"), sep = "\t", header=TRUE)
iddiff <- diff(pg$id)  # same as : pg$id[2:nrow(pg)] - pg$id[1:(nrow(pg)-1)]
iddiff |> table()


download.file("https://gitlab.bfabric.org/wolski/prolfquadata/-/raw/master/inst/quantdata/MSFragger_IonStar2018_PXD003881.zip?inline=false", "MSFragger_IonStar2018_PXD003881.zip")
tmp <- unzip("MSFragger_IonStar2018_PXD003881.zip", list = TRUE)
