Class 7: R functions and packages
================
Joan M. Valls Cuevas
4/24/2019

More on function writting
-------------------------

First we will revisit out function from the last lecture

``` r
source("http://tinyurl.com/rescale-R")
```

Test the **rescale()** function

``` r
rescale(1:10)
```

    ##  [1] 0.0000000 0.1111111 0.2222222 0.3333333 0.4444444 0.5555556 0.6666667
    ##  [8] 0.7777778 0.8888889 1.0000000

Test **rescale()** using a string

``` r
##rescale(c(1:10, "string"))
```

Now that we see an error was thrown to us, we will use a **warning()** or **stop()** function. use the rescale2() function which applies this

``` r
##rescale2(c(1:10, "string"))
```

Function practice
-----------------

Write a function to identify NA elements in two vectors

Starting with a simple example to test body of function

``` r
x <- c(1,2,NA,3,NA)
y <- c(NA, 3, NA, 3, 4)
is.na(x)
```

    ## [1] FALSE FALSE  TRUE FALSE  TRUE

``` r
is.na(x)
```

    ## [1] FALSE FALSE  TRUE FALSE  TRUE

Test the logical AND operator to see if it detects both NA values

``` r
is.na(x) & is.na(y)
```

    ## [1] FALSE FALSE  TRUE FALSE FALSE

Find number of TRUE statements that are both TRUE

``` r
sum(is.na(x) & is.na(y))
```

    ## [1] 1

Using code snippet to make function

``` r
both_na <- function(x,y) {
  sum(is.na(x) & is.na(y))
} 
```

apply function that is made

``` r
both_na(x,y)
```

    ## [1] 1

second test case

``` r
both_na(c(NA,NA,NA), c(1,NA,NA,NA))
```

    ## Warning in is.na(x) & is.na(y): longer object length is not a multiple of
    ## shorter object length

    ## [1] 3

Need to check if length is the same of both vector.

``` r
x <- c(NA,NA,NA)
y <- c(1,NA,NA,NA)
length(x) != length(y)
```

    ## [1] TRUE

Try the both\_na3() function with extra features

``` r
both_na3 <- function(x, y) {

if(length(x) != length(y)) {
stop("Input x and y should be vectors of the same length")

}
  na.in.both <- ( is.na(x) & is.na(y) )
  na.number  <- sum(na.in.both)
  na.which   <- which(na.in.both)
  message("Found ", na.number, " NA's at position(s):",
          paste(na.which, collapse=", ") )
  return( list(number=na.number, which=na.which) )
}
```

Write function to grade homework
================================

start with a simple example

``` r
#student 1
stu1 <- c(100, 100, 100, 100, 100, NA ,100, 90)
#student 2
stu2 <- c(100,NA,90,90,90,90,87,80)

min(stu1)
```

    ## [1] NA

``` r
which.min(stu1)
```

    ## [1] 8

``` r
NA_val <- NA
stu1[which.min(stu1)] <- NA_val
stu1
```

    ## [1] 100 100 100 100 100  NA 100  NA

``` r
grade <- function(x) {

(sum(x, na.rm = TRUE) - min(x ,  na.rm = TRUE)) / (length(x) - 1)


}
```

``` r
grade(stu1)
```

    ## [1] 71.42857

``` r
url <- "https://tinyurl.com/gradeinput"

students <- read.csv(url, row.names = 1)

ans <- apply(students, 1, grade)
```

``` r
sort(ans, decreasing = TRUE)
```

    ##  student-7  student-8 student-13  student-1 student-12 student-16 
    ##      94.00      93.75      92.25      91.75      91.75      89.50 
    ##  student-6  student-5 student-17  student-9 student-14 student-11 
    ##      89.00      88.25      88.00      87.75      87.75      86.00 
    ##  student-3 student-19 student-20  student-2 student-18  student-4 
    ##      84.25      82.75      82.75      82.50      72.75      66.00 
    ## student-15 student-10 
    ##      62.50      61.00

one last function example
-------------------------

Find the intersection of two vectors

``` r
input <-  source("http://tinyurl.com/rescale-R")
```

``` r
x <- df1$IDs
y <- df2$IDs

intersect(x,y)
```

    ## [1] "gene2" "gene3"

``` r
x
```

    ## [1] "gene1" "gene2" "gene3"

``` r
y
```

    ## [1] "gene2" "gene4" "gene3" "gene5"

``` r
x %in% y
```

    ## [1] FALSE  TRUE  TRUE

``` r
x[x %in% y]
```

    ## [1] "gene2" "gene3"

``` r
gene_intersect <- function(x, y) {
  cbind(x[x %in% y],
        y[y %in% x])
  
}
```

``` r
merge(df1, df2, by="IDs")
```

    ##     IDs exp.x exp.y
    ## 1 gene2     1    -2
    ## 2 gene3     1     1
