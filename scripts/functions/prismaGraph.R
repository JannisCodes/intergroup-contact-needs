prismaGr <- function (found, found_other, no_dupes, screened, screen_exclusions, 
          full_text, full_text_exclusions, qualitative, quantitative = NULL, 
          labels = NULL, extra_dupes_box = FALSE, title = "", ..., dpi = 72, font_size = 10) {
  
  dupes <- found + found_other - no_dupes
  labels <- list(found = paste("Records identified through", "database searching", sprintf("(n = %d)", found), sep = "\n"), 
                 found_other = paste("Additional records identified", "through other sources", sprintf("(n = %d)", found_other), sep = "\n"), 
                 no_dupes = paste("Records after duplicates removed", sprintf("(n = %d)", no_dupes), sep = "\n"), 
                 dupes = paste("Duplicates excluded", sprintf("(n = %d)", dupes), sep = "\n"), 
                 screened = paste("Records screened", sprintf("(n = %d)", screened), sep = "\n"), 
                 screen_exclusions = paste("Records excluded", "[reasons, see Table 2]", sprintf("(n = %d)", screen_exclusions), sep = "\n"), 
                 full_text = paste("Full-text articles assessed", "for eligibility", sprintf("(n = %d)", full_text), sep = "\n"), 
                 full_text_exclusions = paste("Full-text articles excluded", "[reasons, see Table 2]", sprintf("(n = %d)", full_text_exclusions), sep = "\n"), 
                 qualitative = paste("Studies included in qualitative synthesis",sprintf("(n = %d)", qualitative), sep = "\n"), 
                 quantitative = paste("Studies included in", "quantitative synthesis", "(empirical)", sprintf("(n = %d)", quantitative), sep = "\n"))
  
  dupes_box <- sprintf("nodups -> incex;\n    nodups [label=\"%s\"];", 
                       labels$no_dupes)
  if (extra_dupes_box) 
    dupes_box <- sprintf("nodups -> {incex; dups};\n       nodups [label=\"%s\"];\n       dups [label=\"%s\"]; {rank=same; nodups dups}", 
                         labels$no_dupes, labels$dupes)
  dot_template <- paste0("digraph prisma {labelloc=\"t\"; label=\"", title, "\n\n\";\n    node [shape=\"box\", fontsize = %d];\n    graph [splines=ortho, nodesep=1, dpi = %d]\n    a -> nodups;\n    b -> nodups;\n    a [label=\"%s\"];\n    b [label=\"%s\"]\n    %s\n    incex -> {ex; ft}\n    incex [label=\"%s\"];\n    ex [label=\"%s\"];\n    {rank=same; incex ex}\n    ft -> {qual; ftex};\n    ft [label=\"%s\"];\n    {rank=same; ft ftex}\n    ftex [label=\"%s\"];\n    qual -> quant\n    qual [label=\"%s\"];\n    quant [label=\"%s\"];\n  }")
  
  DiagrammeR::grViz(sprintf(dot_template, font_size, dpi, labels$found, labels$found_other, 
                            dupes_box, labels$screened, labels$screen_exclusions, 
                            labels$full_text, labels$full_text_exclusions, labels$qualitative, 
                            labels$quantitative))
}

prismaGrScales <- function (found, found_other, no_dupes, 
                            #screened, screen_exclusions, 
                            full_text, full_text_exclusions, qualitative, quantitative = NULL, 
                            labels = NULL, extra_dupes_box = FALSE, title = "", ..., dpi = 72, font_size = 10) {
  
  dupes <- found + found_other - no_dupes
  labels <- list(found = paste("Records identified through", "past reviews", sprintf("(n = %d)", found), sep = "\n"), 
                 found_other = paste("Additional records identified", "through own review", sprintf("(n = %d)", found_other), sep = "\n"), 
                 no_dupes = paste("Records after duplicates removed", sprintf("(n = %d)", no_dupes), sep = "\n"), 
                 dupes = paste("Duplicates excluded", sprintf("(n = %d)", dupes), sep = "\n"), 
                 #screened = paste("Records screened", sprintf("(n = %d)", screened), sep = "\n"), 
                 #screen_exclusions = paste("Records excluded", "[reasons, see Table 2]", sprintf("(n = %d)", screen_exclusions), sep = "\n"), 
                 full_text = paste("Full-text articles assessed", "for eligibility", sprintf("(n = %d)", full_text), sep = "\n"), 
                 full_text_exclusions = paste("Full-text articles excluded", "[reasons, see Table 2]", sprintf("(n = %d)", full_text_exclusions), sep = "\n"), 
                 qualitative = paste("Scales included in qualitative synthesis",sprintf("(n = %d)", qualitative), sep = "\n"), 
                 quantitative = paste("Scales included in", "quantitative synthesis", "(empirical)", sprintf("(n = %d)", quantitative), sep = "\n"))
  
  dupes_box <- sprintf("nodups -> ft;\n    nodups [label=\"%s\"];", 
                       labels$no_dupes)
  if (extra_dupes_box) 
    dupes_box <- sprintf("nodups -> {ft; dups};\n       nodups [label=\"%s\"];\n       dups [label=\"%s\"]; {rank=same; nodups dups}", 
                         labels$no_dupes, labels$dupes)
  dot_template <-  paste0("digraph prisma {labelloc=\"t\"; label=\"", title, "\n\n\";\n    node [shape=\"box\", fontsize = %d];\n    graph [splines=ortho, nodesep=1, dpi = %d]\n    a -> nodups;\n    b -> nodups;\n    a [label=\"%s\"];\n    b [label=\"%s\"]\n    %s\n    ft -> {qual; ftex};\n    ft [label=\"%s\"];\n    {rank=same; ft ftex}\n    ftex [label=\"%s\"];\n    qual -> quant\n    qual [label=\"%s\"];\n    quant [label=\"%s\"];\n  }")
  
  DiagrammeR::grViz(sprintf(dot_template, font_size, dpi, labels$found, labels$found_other, 
                            dupes_box, 
                            labels$full_text, labels$full_text_exclusions, labels$qualitative, 
                            labels$quantitative))
}