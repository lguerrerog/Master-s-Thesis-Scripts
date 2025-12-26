# Read the file (no header, separated by whitespace)
cv <- read.table("/Users/macbook/Documents/Master/Master_Thesis/Results/ADMIXTURE/CV_error2.txt", header = FALSE)

# Inspect
head(cv)
#   V1 V2     V3   # V1 = "CV", V2 = K value, V3 = error

# Scatterplot
plot(cv$V2, cv$V3,
     xlab = "K",
     ylab = "CV error",
     main = "Cross-Validation Error vs K",
     pch = 19)      # filled circles
lines(cv$V2, cv$V3, col = "blue") 

# Add labels only where K = 7 or 9
idx <- cv$V2 %in% c(7, 9)
text(cv$V2[idx], cv$V3[idx],
     labels = round(cv$V3[idx], 4),
     pos = 4, cex = 0.55)
