
library(ggplot2)
library(ggrepel)
library(dplyr)
library(stringr)
library(scales)

#set working directory 
setwd("/Users/macbook/Documents/Master/Master_Thesis/Results/PCA")

fn <- "/Users/macbook/Documents/Master/Master_Thesis/Results/PCA/coordinates_1240Kv62__C_sgl_sgh.evec" 
evecDat <- read.table(fn, col.names=c("Sample", "PC1", "PC2", "PC3", "PC4", "PC5","PC6", "PC7", "PC8", "PC9", "PC10", "Pop")) 
popGroups <- read.table("/Users/macbook/Documents/Master/Master_Thesis/Results/PCA/name_list_full_sgHigh_ploty.txt", 
                        col.names=c("Pop", "PopGroup")) 
mergedEvecDat = merge(evecDat, popGroups, by="Pop") 

# 6) Flag Umingmak points for highlighting
mergedEvecDat$IsUmingmak <- grepl("^Umingmak", mergedEvecDat$PopGroup)

#7) subset dataset excluding 
exclude_groups <- c("Umingmak_SGB1", "Umingmak_SgnB")
subsetDat <- subset(mergedEvecDat, !(PopGroup %in% exclude_groups))

# ---------------------------
# 1) custom palette
# ---------------------------
my_colors <- c(
  "Cau"   = "#AB6AD5", "Eur"  = "#1F77B4",     "Bering" = "#AEC7E8",   "WSib" = "#ffff00",
  "SEA"   = "#FFBB78", "Oce"  = "#2CA02C",     "ME"     = "#98DF8A",   "Afr"  = "#D62728",
  "CA"    = "#FF9896", "Ind"  = "#8e3a59",     "SSib"   = "#8C564B",   "SAM"  = "#C49C94",
  "Ath"   = "#E377C2", "NAM"  = "#F7B6D2",     "EA"     = "#7F7F7F",   "ESib" = "#C7C7C7",
  "Russia_Eur"       = "#BCBD22",
  "Umingmak_CnB"     = "#DBDB8D",
  "Umingmak_CB"      = "#17BECF",
  "Umingmak_SgnB"    = "#9EDAE5",
  "Umingmak_SGB1"    = "#DBA13A",
  "Umingmak_SgB"    = "#FF2700"
)

# Ensure PopGroup text is clean and matches color names exactly
mergedEvecDat <- mergedEvecDat |>
  mutate(PopGroup = str_trim(PopGroup))

# Warn if any PopGroup has no color defined (helps avoid blank plots)
missing_cols <- setdiff(unique(mergedEvecDat$PopGroup), names(my_colors))
if (length(missing_cols) > 0) {
  message("These PopGroup values have no color in my_colors: ", paste(missing_cols, collapse=", "))
}

# ---------------------------
# 2) Highlight Umingmak points
# ---------------------------

# label only the highlighted points
label_um <- TRUE
lab_data <- subsetDat |> filter(IsUmingmak)

# ---------------------------
# 3) Load evals file and set axis names
# ---------------------------
pc1_lab <- "PC1"
pc2_lab <- "PC2"

# Load evals file:
 evals <- readr::read_table("/Users/macbook/Documents/Master/Master_Thesis/Results/PCA/variance_1240Kv62__C_sgl_sgh.eval",
                            col_names = "eig")
 
 # ---------------------------
 # 4) Build the figure
 # ---------------------------
 base_size <- 12  
 p <- ggplot(subsetDat, aes(PC1, PC2)) +
   # zero lines
   geom_hline(yintercept = 0, linewidth = 0.25, linetype = "dashed", color = "grey60") +
   geom_vline(xintercept = 0, linewidth = 0.25, linetype = "dashed", color = "grey60") +
   
   # main points
   geom_point(aes(color = PopGroup), size = 1.9, alpha = 0.85, stroke = 0) +
   
   # highlight layer (larger, black outline)
   geom_point(
     data = subset(subsetDat, IsUmingmak),
     aes(color = PopGroup),
     size = 4, 
     shape = 21, 
     fill = NA,         # keep transparent center
     stroke = 1.2       # width of the outline
   ) +
   
   scale_color_manual(values = my_colors, guide = guide_legend(override.aes = list(size = 3))) +
   
   labs(
     x = pc1_lab,
     y = pc2_lab,
     color = "Population group",
     #title = "Worldwide PCA",
     #subtitle = "Custom colors with Umingmak samples highlighted",
     caption = NULL
   ) +
   
   coord_equal() +
   theme_minimal(base_size = base_size) +
   theme(
     plot.title      = element_text(face = "bold"),
     plot.subtitle   = element_text(margin = margin(b = 6)),
     legend.position = "right",
     legend.title    = element_text(),
     legend.text     = element_text(size = base_size - 1),
     panel.grid.minor = element_blank(),
     panel.grid.major = element_line(linewidth = 0.25, color = "grey90"),
     axis.title.x    = element_text(margin = margin(t = 8)),
     axis.title.y    = element_text(margin = margin(r = 8))
   )
 
 # Print to screen
 print(p)
 
 # ---------------------------
 # 5) Export
 # ---------------------------
 # Output direction:
 out_dir <- "/Users/macbook/Documents/Master/Master_Thesis/Results/PCA"
 
 # High-res PNG (raster)
 ggsave(file.path(out_dir, "PCA_worldwide_pretty.png"), p,
        width = 9, height = 6, dpi = 600, bg = "white")
 
 # PDF (vector)
 ggsave(file.path(out_dir, "PCA_worldwide_pretty.pdf"), p,
        width = 9, height = 6, device = cairo_pdf, bg = "white")
 
 # SVG (vector)
 ggsave(file.path(out_dir, "PCA_worldwide_pretty.svg"), p,
        width = 9, height = 7, device = svglite::svglite, bg = "white")
 
