#' Read HDF5 files
#' @param file path to the HDF5 file to be read
#' @param name name of the attr
#' @param ... parameters passed to \code{\link[rhdf5]{h5read}}
#' @importFrom rhdf5 h5ls h5read
#' @importFrom dplyr filter_ %>%
#' @return an array of class \code{brain}
#' @export
#' @examples
#' file <- "~/Data/barresi/AT_1_Probabilities.h5"
#' raw <- read_h5(file)
#' class(raw)
#' dim(raw)

read_h5 <- function(file, name = NULL, ...) {
  objs <- rhdf5::h5ls(file) %>%
    filter_(~otype == "H5I_DATASET")
  if (is.null(name)) {
    # take the first HDF5 dataset
    name <- objs$name[1]
  }
  x <- rhdf5::h5read(file, name, ...)
  if (length(dim(x)) > 3) {
    # take the second image
    x <- x[2,,,]
  }
  class(x) <- append("brain", class(x))
  return(x)
}

#' Tidy 3D brain image data
#' @importFrom dplyr %>% mutate_ select_ as.tbl
#' @importFrom broom tidy
#' @export
#' @examples
#' file <- "~/Data/barresi/AT_1_Probabilities.h5"
#' if (require(dplyr)) {
#'   tidy_brain <- file %>%
#'     read_h5() %>%
#'     tidy()
#'   class(tidy_brain)
#'   glimpse(tidy_brain)
#' }

tidy.brain <- function(x, ...) {
  res <- x %>%
    as.data.frame.table() %>%
    mutate_(x = ~as.integer(Var1),
            y = ~as.integer(Var2),
            z = ~as.integer(Var3)) %>%
    select_(~x, ~y, ~z, ~Freq) %>%
    as.tbl()
  class(res) <- append("tbl_brain", class(res))
  return(res)
}

#' Print a slice of a 3D image
#' @export
#' @examples
#' file <- "~/Data/barresi/AT_1_Probabilities.h5"
#' raw <- read_h5(file)
#' image(raw)
#' image(raw, z = 23)


image.brain <- function(x, z = NULL, ...) {
  if (is.null(z)) {
    z <- floor(dim(x)[3] / 2)
  }
  image(x[,,z])
}

#' Plot a 3B image of a brain
#' @inheritParams rgl::plot3d
#' @param threshold value below which points will not be plotted
#' @importFrom rgl plot3d
#' @importFrom dplyr %>% mutate_ filter_
#' @export
#' @examples
#' file <- "~/Data/barresi/AT_1_Probabilities.h5"
#' if (require(dplyr)) {
#'   tidy_brain <- file %>%
#'     read_h5() %>%
#'     tidy()
#'   plot3d(tidy_brain)
#' }

plot3d.tbl_brain <- function(x, threshold = 0.9, ...) {
  xyz <- x %>%
    mutate_(gray_val = ~gray(1 - Freq)) %>%
    filter_(~Freq >= threshold)
  if (nrow(xyz) > 250000) {
    warning(paste("You are about to plot", nrow(xyz),
                  "points. Consider increasing the threshold parameter from", threshold))
  }
  rgl::plot3d(xyz, type = "p", alpha = xyz$Freq,
#         width = 900, height = 800,
         xlab = "x", ylab = "y", zlab = "z",
         size = 2, col = xyz$gray_val, ...)
  plot3d_model(xyz)
}

#' @importFrom mosaic makeFun

plot3d_model <- function(xyz, ...) {
  # quadratic plane
  mod1 <- lm(z ~ x + y + I(x^2), data = xyz)
  bent_plane <- mosaic::makeFun(mod1)
  plot3d(bent_plane, alpha = 0.5, col = "dodgerblue",
         xlim = range(xyz$x), ylim = range(xyz$y), zlim = range(xyz$z),
         add = TRUE)

  # flat plane
  mod2 <- lm(y ~ x + z, data = xyz)
  flat_plane <- mosaic::makeFun(mod2)
  plot3d(flat_plane, alpha = 0.5, col = "dodgerblue",
         xlim = range(xyz$x), ylim = range(xyz$y), zlim = range(xyz$z),
         add = TRUE)

  # their intersection
  a <- coef(mod1)["I(x^2)"]
  b <- coef(mod1)["x"]
  c <- coef(mod1)["y"]
  d <- coef(mod1)["(Intercept)"]
  e <- coef(mod2)["x"]
  f <- coef(mod2)["z"]
  g <- coef(mod2)["(Intercept)"]
  a_prime <- (a*f) / (1 - c*f)
  b_prime <- (e + b*f) / (1 - c*f)
  c_prime <- (d*f + g) / (1 - c*f)
  # make a parametric equation in terms of t
  t <- seq(0, 1000)
  xyz_rope <- cbind(x = t,
                    y = a_prime*t^2 + b_prime*t + c_prime,
                    z = (a + c*a_prime)*t^2 + (b + c*b_prime)*t + c*c_prime + d)
  plot3d(xyz_rope, alpha = 0.5, col = "red", size = 5,
         xlim = range(xyz$x), ylim = range(xyz$y), zlim = range(xyz$z),
         add = TRUE)
}




