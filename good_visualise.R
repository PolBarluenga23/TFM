library(tidyverse)
library(DT)
library(htmlwidgets)
library(htmltools)

tsv_files <- list.files(".", pattern = "\\.tsv$", full.names = TRUE)

html_files <- c()

for(file_path in tsv_files){
  
  blast <- read_tsv(file_path, 
                    col_names = c("qseqid","sseqid","pident","length","mismatch",
                                  "gapopen","qstart","qend","sstart","send",
                                  "evalue","bitscore","qseq","sseq"),
                    show_col_types = FALSE) %>%
    mutate(
      evalue = as.numeric(evalue),
      bitscore = as.numeric(bitscore),
      pident = as.numeric(pident),
      query_virus = "WN1",
      subject_virus = str_extract(sseqid, "^[^|]+"),
      query_acc = str_extract(qseqid, "(?<=\\|)[^|]+(?=\\|)"),
      subject_acc = str_extract(sseqid, "(?<=\\|)[^|]+(?=\\|)"),
      region = str_extract(qseqid, "(?<=\\|)[^|]+(?=\\|Flaviviridae)"),
      mismatch_perc = round((mismatch / length) * 100, 1)
    ) %>%
    arrange(desc(bitscore)) %>%
    slice_head(n = 100)
  
  # ---- PERFECT ALIGNMENT PREVIEW ----
  blast <- blast %>%
    rowwise() %>%
    mutate(
      match_line = paste0(
        ifelse(strsplit(qseq, "")[[1]][1:60] ==
                 strsplit(sseq, "")[[1]][1:60],
               "|", " "),
        collapse = ""
      ),
      alignment_preview = paste0(
        "<pre style='margin:0; font-family:monospace;'>",
        substr(qseq, 1, 60), "\n",
        match_line, "\n",
        substr(sseq, 1, 60),
        "</pre>"
      )
    ) %>%
    ungroup()
  
  dt <- datatable(
    blast %>% select(query_virus, subject_virus, pident, bitscore,
                     evalue, region, alignment_preview),
    filter = "top",
    escape = FALSE,
    options = list(
      pageLength = 12,
      order = list(list(3, "desc")),
      searchHighlight = TRUE,
      scrollX = TRUE,
      columnDefs = list(list(width = "260px", targets = 6))
    ),
    rownames = FALSE
  ) %>%
    formatRound(columns = c("pident", "bitscore", "evalue"), digits = 2) %>%
    formatStyle("bitscore", 
                backgroundColor = styleInterval(c(400, 450), 
                                                c("white", "lightyellow", "lightgreen"))) %>%
    formatStyle("pident", 
                backgroundColor = JS("function(row){return row[2]>95?'lightgreen':row[2]>85?'lightyellow':'white';}"))
  
 
  name <- paste0(tools::file_path_sans_ext(basename(file_path)), ".html")
  saveWidget(dt, name, selfcontained = TRUE)
  
  html_files <- c(html_files, name)
}



tabs <- lapply(seq_along(html_files), function(i){
  tags$div(
    class = if(i==1) "tab-pane active" else "tab-pane",
    id = paste0("tab", i),
    tags$iframe(
      src = html_files[i],
      style = "width:100%; height:900px; border:none;"
    )
  )
})

nav <- tags$ul(
  class = "nav nav-tabs",
  lapply(seq_along(html_files), function(i){
    tags$li(
      class = if(i==1) "active" else "",
      tags$a(
        href = paste0("#tab", i),
        `data-toggle` = "tab",
        tools::file_path_sans_ext(html_files[i])
      )
    )
  })
)

page <- tagList(
  tags$head(
    tags$link(rel="stylesheet",
              href="https://maxcdn.bootstrapcdn.com/bootstrap/3.4.1/css/bootstrap.min.css"),
    tags$script(src="https://ajax.googleapis.com/ajax/libs/jquery/3.5.1/jquery.min.js"),
    tags$script(src="https://maxcdn.bootstrapcdn.com/bootstrap/3.4.1/js/bootstrap.min.js")
  ),
  tags$h2("BLAST Alignment"),
  nav,
  tags$div(class="tab-content", tabs)
)

save_html(page, "ALL_WN1_BLAST_TABS.html")
browseURL("ALL_WN1_BLAST_TABS.html")
