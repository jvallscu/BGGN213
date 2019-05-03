#' ---
#' title: "Class 5: R Graphics "
#' author: "Joan M. Valls Cuevas"
#' date: "April 19, 2019"
#' ---


# Class 5 R graphics and plots

# get the data in
weight <- read.table("bimm143_05_rstats/weight_chart.txt", 
                     header = TRUE)

# plot a  scaterplot of age vs weight

plot(weight[,1], weight[,2],
     xlab = "Age (months)",
     ylab = "Weight (Kg)", main = "Age vs. weight", type = "b", pch = 15, cex = 1.3)


# Bar Plot section 2B
feature <- read.table("bimm143_05_rstats/feature_counts.txt",
                      header = TRUE, sep = "\t")


# old par values

old.par <- par()$mar

par(mar=c(5,11,4,3))

barplot(feature$Count, horiz = T, ylab = "A title",
        names.arg = feature$Feature, main = "Features",
        las = 1)
# 3A

mf <- read.delim("bimm143_05_rstats/male_female_counts.txt", header = TRUE )


barplot(mf$Count, names.arg = mf$Sample, las = 2, col = rainbow(10))

rainbow(10)

# 3B

genes <- read.delim("bimm143_05_rstats/up_down_expression.txt")


table(genes$State)
length(genes$State)

palette(c("blue", "grey", "red"))
plot(genes$Condition1, genes$Condition2, col = genes$State, xlab = "Expression condition 1", ylab = "Expression condition 2")



