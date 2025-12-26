##
## Base-R script to plot ADMIXTURE Q matrix
##

# ---------- input ----------
ind_file <- "/Users/macbook/Documents/Master/Master_Thesis/Results/ADMIXTURE/subset2_ADMIXTURE.ind"        # .ind file
q_file   <- "/Users/macbook/Documents/Master/Master_Thesis/Results/ADMIXTURE/subset2_ADMIXTURE.pruned.7.Q" # .Q file
K        <- 7                                    # number of clusters
pdf_out  <- "/Users/macbook/Documents/Master/Master_Thesis/Results/ADMIXTURE/admixtureK7.pdf"      # output figure

my_cols <- c(
  "#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#9BCD9B",
  "#e6ab02", "#4F94CD", "#EE0000", "#8dd3c7", "#fb8072",
  "#CD96CD", "#80b1d3", "#8B5A2B", "#8B0000", "#FFEC8B"
)

# ---------------------------------

# Read the .ind (individual.id, sex, group.id)
ind <- read.table(ind_file, stringsAsFactors = FALSE)
colnames(ind) <- c("id","sex","group")

# Read the Q matrix
Q <- as.matrix(read.table(q_file))

if (nrow(Q) != nrow(ind)) stop("Individuals in .ind and .Q do not match!")

#----------Order of populations--------------

# 1. Read in the desired order
pop_order <- scan("/Users/macbook/Documents/Master/Master_Thesis/Results/ADMIXTURE/reversed_pop_list_ADMIXTURE2.txt", what = character(), quiet = TRUE)

# 2. Make sure every group in ind is in that list
if (!all(unique(ind$group) %in% pop_order)) {
  warning("Some population names in your data are not in pop_order.txt")
}

# 3. Create a factor with the specified levels
ind$group <- factor(ind$group, levels = pop_order)

# 4. Order individuals by that factor
ord <- order(ind$group)
ind <- ind[ord, ]
Q   <- Q[ord, ]

#-----------Plotting---------------
# Colours
if (length(my_cols) < K) {
  my_cols <- rep(my_cols, length.out = K)
}

# Open PDF
pdf(pdf_out, width = 7, height = max(6, nrow(Q) * 0.15))
oldpar <- par(no.readonly = TRUE)     # save current settings
par(mar = c(6, 12, 2, 2))             # bottom, left, top, right


# Set up the plot
bar_mid <- barplot(
  t(Q),                          # transpose so clusters stack
  col = my_cols,
  border = NA,
  space = 0.02,
  horiz = TRUE,
  xlab = "Ancestry proportion",
  ylab = "",
  xaxt = "s", yaxt = "n",                      # x-axis labels manually added
  xlim   = c(0, 1.15)   # allow extra space
)

# 1. group sizes in the plotting order. Add population group separators and labels
pop_sizes <- table(ind$group)
#pop_names <- names(pop_sizes)

# 2. find the midpoints returned by barplot for each individual
cum_sizes <- cumsum(pop_sizes)

# 3. get the coordinate between the last bar of one group and the first of the next
#   halfway point between bar_mid[i] and bar_mid[i+1]
boundary_pos <- (bar_mid[cum_sizes[-length(cum_sizes)]] +
                   bar_mid[cum_sizes[-length(cum_sizes)] + 1]) / 2

# 4. draw horizontal separator lines
# draw the horizontal lines for populations
for (pos in boundary_pos) {
  segments(x0 = 0, x1 = 1,    # draw only from x=0 to x=1
           y0 = pos, y1 = pos,
           col = "#333333",
           lwd = 0.03)
}

# Add population labels in the middle
# compute midpoint for each group
group_centers <- sapply(seq_along(pop_sizes), function(i) {
  idx_start <- if (i == 1) 1 else sum(pop_sizes[1:(i - 1)]) + 1
  idx_end   <- sum(pop_sizes[1:i])
  mean(bar_mid[idx_start:idx_end])            # average of individual midpoints
})
axis(2, at = group_centers, labels = names(pop_sizes), las = 1, cex.axis = 0.6, col = "#333333")
#----------Other optioins of labels-----------
#Add population labels at the first position
# position of *first* individual in each group
#first_pos <- c(1, head(cum_sizes + 1, -1))
#axis(2, at = bar_mid[first_pos], labels = names(pop_sizes), las = 1, cex.axis = 0.7)

#Adds population name on every bar
#axis(2,
#  at     = bar_mid,            # same length as number of individuals
#  labels = ind$group,          # group name for each individual
#  las    = 1,                  # horizontal text
#  cex.axis = 0.6               # shrink if you have many samples
#)

#----------Geographical group axis-----------
## 1. Read the mapping file
##    File format: pop<TAB>geo_group
geo_map <- read.table("/Users/macbook/Documents/Master/Master_Thesis/Results/ADMIXTURE/pop_list_geo_ADMIXTURE.txt",
                      header = FALSE,
                      sep = "\t",
                      stringsAsFactors = FALSE,
                      col.names = c("pop", "geo"))

## Ensure populations match
if (!all(names(pop_sizes) %in% geo_map$pop)) {
  warning("Some populations are missing in the geo map!")
}

## 2. Merge to attach the geographic group for each population
pop_df <- data.frame(pop = names(pop_sizes),
                     size = as.vector(pop_sizes),
                     stringsAsFactors = FALSE)

pop_df <- merge(pop_df, geo_map, by = "pop")

## Keep the same plotting order
# pop_df contains: pop, size, geo   (and is in plotting order)
pop_df <- pop_df[match(names(pop_sizes), pop_df$pop), ]



#---------------Final code for geo-labels-----------
# 1. Get index range (start & end) of each geo group

geo_ranges <- tapply(seq_len(sum(pop_df$size)),
                     rep(pop_df$geo, pop_df$size),
                     range)

# 2. Horizontal position for the bracket (outside bars)
x_bracket <- 1.01   

# 3. Draw brackets with labels
for (g in names(geo_ranges)) {
  idx <- geo_ranges[[g]]
  y1  <- bar_mid[min(idx)]   # top of group
  y2  <- bar_mid[max(idx)]   # bottom of group
  
  # vertical main line of the bracket
  segments(x_bracket, y1, x_bracket, y2, lwd = 2, col="#4D4D4D")
  
  # small horizontal ticks at top and bottom
  tick_len <- 0.01           # length of the horizontal tick
  segments(x_bracket, y1, x_bracket + tick_len, y1, lwd = 2, col = "#4D4D4D")
  segments(x_bracket, y2, x_bracket + tick_len, y2, lwd = 2, col = "#4D4D4D")
  
  # vertical text label centered on the bracket
  text(x_bracket + tick_len + 0.02,
       mean(c(y1, y2)),
       labels = g,
       srt = 270,             # rotate text
       adj = 0.5,
       cex = 0.6)
}

dev.off()

cat("Saved plot to", pdf_out, "\n")
