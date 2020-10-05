
# weitrix 1.1.8

* Use geom_bin2d in weitrix_calplot scatterplots.


# weitrix 1.1.7

* Add mu_min, mu_max arguments to weitrix_calibrate_all.


# weitrix 1.1.6

* Use glm2, which is less prone to optimization failure.
* Auto-disable parallel processing if X11 device is open.
* Add SLAM-Seq vignette.


# weitrix 1.1.5

* well_knotted_spline for natural splines with good choice of knots.


# weitrix 1.1.4

* weitrix_confects for differential testing.
* weitrix_rms_confects now called weitrix_sd_confects.
* weitrix_calplot now shows mean trend and mean +/- standard deviation trend.


# weitrix 1.1.3

* weitrix_calibrate_all now includes a simple scaling factor to account for
  residuals from a fitted model being smaller than residuals from the true 
  model.
* weitrix_calplot blue guidelines are similarly adjusted.
* Add weitrix_rms_confects to find rows with confidently excessive variation.


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
