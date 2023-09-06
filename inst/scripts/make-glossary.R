# This script converts the glossary stored as an excel file in
# "inst/extdata/slo_glo_table.xlsx" into a data frame saved in the "data" folder

# Load packages
library(readxl)


# Read glossary
glossary <- read_excel("inst/extdata/slo_glo_table.xlsx")
glossary <- as.data.frame(glossary)


# Save
save(glossary, file =  "data/glossary.rda")
