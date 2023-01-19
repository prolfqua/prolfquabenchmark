
tex <- readLines("paper.tex")

abstract <- tex[ grep("begin\\{abstract\\}", tex):(grep("end\\{abstract\\}", tex))]
writeLines(abstract, "../../../prolfqua-achems/abstract.tex")

introduction <- tex[ grep("hypertarget\\{introduction\\}", tex):(grep("hypertarget\\{methods\\}", tex) - 1)]
writeLines(introduction, "../../../prolfqua-achems/introduction.tex")

methods <- tex[ grep("hypertarget\\{methods\\}", tex):(grep("hypertarget\\{results-and-discussion\\}", tex) - 1)]
writeLines(methods, "../../../prolfqua-achems/methods.tex")

# create results.tex
results <- tex[ grep("hypertarget\\{results-and-discussion\\}", tex):(grep("hypertarget\\{conclusion\\}", tex) - 1)]

sB <- grep("begin\\{Shaded\\}", results)
sE <- grep("end\\{Shaded\\}", results)

results <- results[-(sB[1]:sE[2])]
results <- append(results, "\\TODO{codesnippets/codesnippet2.pdf}", after = (sB[1] - 1))

sB <- grep("begin\\{Shaded\\}", results)
sE <- grep("end\\{Shaded\\}", results)
results <- results[-(sB[1]:sE[1])]
results <- append(results, "\\TODO{codesnippets/codesnippet3.pdf}", after = (sB[1] - 1))

writeLines(results, "../../../prolfqua-achems/results.tex")

conclusion <- tex[ grep("hypertarget\\{conclusion\\}", tex):(grep("hypertarget\\{acknowledgements\\}", tex) - 1)]
writeLines(conclusion, "../../../prolfqua-achems/conclusion.tex")

abbreviations <- tex[ grep("hypertarget\\{abbreviations\\}", tex):(grep("hypertarget\\{supporting-information\\}", tex) - 1)]
writeLines(abbreviations, "../../../prolfqua-achems/abbreviations.tex")

abbreviations <- tex[ grep("hypertarget\\{supporting-information\\}", tex):(grep("hypertarget\\{references\\}", tex) - 1)]
writeLines(abbreviations, "../../../prolfqua-achems/supporting.tex")
