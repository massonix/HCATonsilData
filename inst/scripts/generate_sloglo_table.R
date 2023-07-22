library("googlesheets4")

gsheet_url <- "https://docs.google.com/spreadsheets/d/1myzlfmyh8QjSha3-c3Z32QtoKUhffdjx0A3rj1wJzac/edit#gid=0"

sloglo_df <- read_sheet(gsheet_url) |>
  as.data.frame()
rownames(sloglo_df) <- sloglo_df$annotation_detailed

sloglo_df

DT::datatable(sloglo_df)

write.table(sloglo_df,
  file = "inst/extdata/sloglo_tabular.csv",
  col.names = TRUE,
  row.names = TRUE
)

# check it reads in correctly
read.table("inst/extdata/sloglo_tabular.csv")
