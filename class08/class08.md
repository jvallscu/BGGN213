Class08
================
Joan M. Valls Cuevas
4/26/2019

K means clustering and unsupervised learning
--------------------------------------------

data split up into a pre-defined number of clusters

lots of calculations.

**kmeans(x, center = (num of groups), nstart = (num of iterations))**

scree plots, plot

Some exampples for clustering
=============================

``` r
#Generate some example data for clustering 
tmp <- c(rnorm(30,-3), rnorm(30,3))
x <- cbind(x = tmp, y = rev(tmp))

km <- kmeans(x, centers = 2, nstart = 20)

plot(x, col = km$cluster)
```

![](class08_files/figure-markdown_github/unnamed-chunk-1-1.png)

``` r
table(km$cluster)
```

    ## 
    ##  1  2 
    ## 30 30

plot x colored by the kmeans cluster assignment

``` r
plot(x, col = km$cluster)
points(km$centers, col = "blue", pch = 18, cex = 8)
```

![](class08_files/figure-markdown_github/unnamed-chunk-2-1.png)

Hierarchical clustering
=======================

Here we dont have to spell out the number of clusters before hand but we have to give a it a distance matrix

bottom up or top/down types

``` r
d <- dist(x)
hc <- hclust(d)
```

Let's plot the resutls

``` r
plot(hc)
abline( h = 6, col = "red")
```

![](class08_files/figure-markdown_github/unnamed-chunk-4-1.png)

``` r
cutree(hc, h = 6)
```

    ##  [1] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2
    ## [36] 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2

``` r
gp2 <- cutree(hc, k=2)
gp3 <- cutree(hc, k=3)
table(gp2)
```

    ## gp2
    ##  1  2 
    ## 30 30

``` r
table(gp3)
```

    ## gp3
    ##  1  2  3 
    ## 30  8 22

``` r
# Step 1. Generate some example data for clustering
x2 <- rbind(
  matrix(rnorm(100, mean=0, sd = 0.3), ncol = 2),   # c1
  matrix(rnorm(100, mean = 1, sd = 0.3), ncol = 2), # c2
  matrix(c(rnorm(50, mean = 1, sd = 0.3),           # c3
           rnorm(50, mean = 0, sd = 0.3)), ncol = 2))
colnames(x2) <- c("x", "y")
# Step 2. Plot the data without clustering
plot(x2)
```

![](class08_files/figure-markdown_github/unnamed-chunk-5-1.png)

``` r
# Step 3. Generate colors for known clusters
#         (just so we can compare to hclust results)
col <- as.factor( rep(c("c1","c2","c3"), each=50) )
plot(x2, col=col)
```

![](class08_files/figure-markdown_github/unnamed-chunk-5-2.png)

Q. Use the dist(), hclust(), plot() and cutree() functions to return 2 and 3 clusters Q. How does this compare to your known 'col' groups?

``` r
dist_c <- dist(x2)
hc2 <- hclust(dist_c)
plot(hc2)
abline(h = 2, col = "red")
```

![](class08_files/figure-markdown_github/unnamed-chunk-6-1.png)

``` r
c1 <- cutree(hc2, k = 2)
c2 <- cutree(hc2, k = 3)
table(c1)
```

    ## c1
    ##   1   2 
    ##  48 102

``` r
table(c2)
```

    ## c2
    ##  1  2  3 
    ## 48 48 54

``` r
plot(x2, col = c1)
```

![](class08_files/figure-markdown_github/unnamed-chunk-6-2.png)

``` r
plot(x2, col = c2)
```

![](class08_files/figure-markdown_github/unnamed-chunk-6-3.png) \#\# Principal component Analysis (PCA)

We will use the base R **prcom()** function

``` r
## You can also download this file from the class website!
mydata <- read.csv("https://tinyurl.com/expression-CSV",
row.names=1)
head(mydata)
```

    ##        wt1 wt2  wt3  wt4 wt5 ko1 ko2 ko3 ko4 ko5
    ## gene1  439 458  408  429 420  90  88  86  90  93
    ## gene2  219 200  204  210 187 427 423 434 433 426
    ## gene3 1006 989 1030 1017 973 252 237 238 226 210
    ## gene4  783 792  829  856 760 849 856 835 885 894
    ## gene5  181 249  204  244 225 277 305 272 270 279
    ## gene6  460 502  491  491 493 612 594 577 618 638

There are 100 genes in this dataset

``` r
pca <- prcomp( t(mydata), scale = TRUE)
summary(prcomp( t(mydata), scale = TRUE))
```

    ## Importance of components:
    ##                           PC1    PC2     PC3     PC4     PC5     PC6
    ## Standard deviation     9.6237 1.5198 1.05787 1.05203 0.88062 0.82545
    ## Proportion of Variance 0.9262 0.0231 0.01119 0.01107 0.00775 0.00681
    ## Cumulative Proportion  0.9262 0.9493 0.96045 0.97152 0.97928 0.98609
    ##                            PC7     PC8     PC9      PC10
    ## Standard deviation     0.80111 0.62065 0.60342 3.348e-15
    ## Proportion of Variance 0.00642 0.00385 0.00364 0.000e+00
    ## Cumulative Proportion  0.99251 0.99636 1.00000 1.000e+00

``` r
attributes(pca)
```

    ## $names
    ## [1] "sdev"     "rotation" "center"   "scale"    "x"       
    ## 
    ## $class
    ## [1] "prcomp"

``` r
plot(pca$x[,1], pca$x[,2])
```

![](class08_files/figure-markdown_github/unnamed-chunk-10-1.png)

``` r
pca.var <- pca$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100,1)
pca.var.per
```

    ##  [1] 92.6  2.3  1.1  1.1  0.8  0.7  0.6  0.4  0.4  0.0

``` r
xlab <- paste("PC1 (", pca.var.per[1],"%)", sep="")
ylab <- paste("PC1 (", pca.var.per[2],"%)", sep="")
```

``` r
plot(pca$x[,1],pca$x[,2], xlab=xlab, ylab=ylab)
text(pca$x[,1],pca$x[,2], colnames(mydata))
```

![](class08_files/figure-markdown_github/unnamed-chunk-13-1.png)

PCA of UK data sets
===================

``` r
UK_data <- read.csv("UK_foods.csv", row.names = 1)

dim(UK_data)
```

    ## [1] 17  4

``` r
#View(UK_data)
barplot(as.matrix(UK_data), beside = T, col=rainbow(nrow(UK_data)))
```

![](class08_files/figure-markdown_github/unnamed-chunk-15-1.png)

``` r
pairs(UK_data, col=rainbow(10), pch=16)
```

![](class08_files/figure-markdown_github/unnamed-chunk-16-1.png)

``` r
# Use the prcomp() PCA function 
pca <- prcomp( t(UK_data) )
summary(pca)
```

    ## Importance of components:
    ##                             PC1      PC2      PC3       PC4
    ## Standard deviation     324.1502 212.7478 73.87622 4.189e-14
    ## Proportion of Variance   0.6744   0.2905  0.03503 0.000e+00
    ## Cumulative Proportion    0.6744   0.9650  1.00000 1.000e+00

``` r
plot(pca$x[,1], pca$x[,2])
text(pca$x[,1], pca$x[,2], colnames(UK_data))
```

![](class08_files/figure-markdown_github/unnamed-chunk-18-1.png)
