
library(weitrix)
library(usethis)

set.seed(1)

scores <- c(1,2,3,4,5)
loadings <- c(-1,0,1,0,-1,0,1)

n <- 100
rows <- sample.int(n, n=length(loadings), replace=TRUE)
cols <- sample.int(n, n=length(scores), replace=TRUE)
vals <- scores[cols]*loadings[rows] + rnorm(n)  

x <- tapply(vals, list(rows,cols), mean)
w <- tapply(vals, list(rows,cols), length, default=0)

simwei <- as_weitrix(x, w)
colData(simwei)$true_score <- scores

use_data(simwei, overwrite=TRUE)

