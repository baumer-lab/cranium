An R package to quantify radial bridge images

Installation
------------

`cranium` depends on the`hdf5r` package or `rhdf5` package for reading HDF5 files. You could download `hdf5r` from CRAN. However, `rhdf5` package is not available on CRAN--it is a bioconductor package. You must install it from Bioconductor. To do this, you have to run the following code once:

``` r
install.packages("BiocInstaller",
                 repos = "http://bioconductor.org/packages/3.4/bioc")
source("https://bioconductor.org/biocLite.R")
biocLite("rhdf5")
```

If you are having trouble downloading `rhdf5` with the code above, an alternative option would be downloading `rhdf5` from bioconductor using the `pak` package: https://github.com/r-lib/pak

After installing `rhdf5`, you can install `cranium` from GitHub using the following commands:

``` r
install.packages("devtools")
devtools::install_github("beanumber/cranium")
```

Plot a 3D image of a brain
--------------------------

Now, you can simply point to a raw HDF5 file, and quickly render a 3D image of the data, along with our model of it.

``` r
file <- "~/Data/barresi/AT_1_Probabilities.h5"
library(cranium)
library(tidyverse)
tidy_brain <- file %>%
  read_h5() %>%
  tidy()
plot3d(tidy_brain)
```
