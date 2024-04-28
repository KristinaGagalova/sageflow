#' Plot insert loci
#'
#' The function uses a gff file, assigns instert group loci and plots them based on coordinates
#' @param gff gz file as output from sageflow
#' @param position_ref (int) - position of th reference sequence to use for comparison. Default=1
#' @import dplyr, ggplot2, gggenes, bedtoolsr, Biostrings, sys, paletteer, fs, tidyverse, glue
#' @keywords insert loci
#' @export
#' @return plot object ready to be plot
#' @examples
#' plotInsert("/path/to/EVA035_int.fasta.gff3.gz")

assignGroups <- function(gzip.gff,gzip_f=T,id_th=0,cov_th=0,DIST=10000,strain = '',my_path=''){
  #-----------------------------------
  #const bedtools
  #DIST = 10000 #set this value to the overlap (upstream and downstream) to assign groups for plotting
  #const filt - provided as default parameters in the function argument
  #ID = 99
  #COV = 30
  #----------------------------------
  packages <- c("ggplot2","dplyr","bedtoolsr","Biostrings","sys","glue","fs","tidyverse")

  # Install packages not yet installed
  installed_packages <- packages %in% rownames(installed.packages())
  if (any(installed_packages == FALSE)) {
   install.packages(packages[!installed_packages])
  }

  # Packages loading
  invisible(lapply(packages, library, character.only = TRUE))
  
  #library(bedtoolsr)
  #if (my_path == ''){
  #options(bedtools.path = "/home/kgagalova/.linuxbrew/bin")
  #} else {
  #options(bedtools.path = my_path)
  #}

  #testthat::test_package("bedtoolsr")
  system("which bedtools")

  if (gzip_f) {
  gff = read.delim2(gzfile(gzip.gff,'rt'),header = F,comment.char = '#')
  } else {
  gff = read.delim2(gzip.gff,header = F,comment.char = '#')
  }

  if (strain != '') {
  gff0 <- gff %>% filter(grepl(strain,V1))
  } else {
  gff0 <- gff
  }

  gff2 = distinct(gff0, V1,V4,V5,V7,.keep_all=T)  
  #remove duplicated exact hits
  tmp_coor = bedtoolsr::bt.merge(i=gff2,d=DIST)
  tmp_coor$num = row.names(tmp_coor)
  int_group = bedtoolsr::bt.intersect(a=gff2,b=tmp_coor,wao=TRUE)
  
  grouped_hits = int_group %>% 
    dplyr::mutate(V9 = str_split(V9, ";")) %>% 
    unnest(cols = c(V9)) %>% 
    separate(V9, c("name", "key"), sep = "=") %>% 
    pivot_wider(names_from = name, values_from = key ) %>% 
    dplyr::mutate(orientation = if_else(V7 == "+", 1, -1),
           group = stringr::str_c(V1, "-", V13),
           id = as.numeric(id),
           cov = as.numeric(cov)) %>%
   dplyr::filter(id > id_th & cov > cov_th) #%>%
  closeAllConnections()
  colnames(grouped_hits)[c(1,4,5)] <- c("scaf","start","stop")
  return(return(grouped_hits[,c("ID","scaf","start","stop","orientation","group","Name","id","cov")]))
}


plotInsert <- function(grouped_hits,insert='',strain=''){
  packages <- c("ggplot2","gggenes","paletteer","ggrepel")
  
  installed_packages <- packages %in% rownames(installed.packages())
  if (any(installed_packages == FALSE)) {
   install.packages(packages[!installed_packages])
  }
  

  if (insert != '') {
  my_groups = grouped_hits %>% filter(Name == insert) %>% pull(group)
    
  grouped_hits2 <- grouped_hits %>% 
	 filter(group %in% my_groups)
  } else {
    grouped_hits2 <- grouped_hits
  }

  if (strain != '') {
  grouped_hits2 <- grouped_hits2[which(grepl(scaf,strain)),]
  }

  toPlot <- grouped_hits2 %>% 
  split(.$group) %>% 
  discard(~nrow(.x) == 1) %>% 
  bind_rows() %>%
  
  ggplot(aes(xmin = start, xmax = stop, y = group, label = Name, forward = orientation, fill = Name)) +
  gggenes::geom_gene_arrow(arrowhead_height = unit(6, "mm"), 
                  arrowhead_width = unit(2, "mm"),
                  arrow_body_height = unit(5, "mm")) +
  gggenes::geom_gene_label( min.size = 4, grow = TRUE, reflow = TRUE, colour = "black",
                   padding.x = unit(0.01, "mm"), padding.y = unit(0.01, "mm")) +
  facet_wrap(~group, scales = "free", ncol = 1) +
  gggenes::theme_genes() + 
  xlab("Coordinates") +
  theme(legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.title.y = element_blank(),
        axis.text = element_text(color="black")) +
  paletteer::scale_fill_paletteer_d("trekcolors::lcars_series") + 
  scale_color_discrete(guide = "none")
  
  return(toPlot)
}
