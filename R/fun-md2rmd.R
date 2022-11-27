
# renv::install("tmcd82070/tex2rmd")
# renv::install("huhuaping/tex2rmd")

require("tex2rmd")
require(stringr)

files_tar <- list.files("mathpix/")
qmd_tar <- str_replace_all(files_tar, "\\.tex","\\.qmd")
#paste0(qmd_tar,collapse = " ")

input_md <- list.files("mathpix/",full.names = T)

i <- 29
for (i in 26:length(input_md)) {
  
  tex2rmd::tex2rmd(infile = input_md[i],
          ext_out = ".qmd",
          dir_img = "images/",
          ext_img = ".jpg",
          head2_only = TRUE, 
          keep_yml =FALSE)
  
  cat(paste0("finished ",i," file convert ", input_md[i] ),
      sep = "\n")
}




