

# weitrix 1.1.2

* Vignettes use weitrix_calibrate_all, demonstrate weitrix_calplot.
* weitrix_calplot now uses sqrt(weight)*residual on y axis.


# weitrix 1.1.1

* Switch calibration from linear model on log dispersions to using a gamma GLM.
* Add weitrix_calibrate_all function for very flexible calibration.
* Add weitrix_calplot to examine quality of weights.


# weitrix 0.99.2

* Automatic package builder wants package to depend on R version 4.0.


# weitrix 0.99.1

* weitrix now uses the BiocParallel default parallel processing, rather than
  the DelayedArray default. (DelayedArray does not use parallel processing
  by default as of version 0.14.0.)


# weitrix 0.99.0

* Bioconductor submission.


# weitrix 0.1.0

* Initial version.
