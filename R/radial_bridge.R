#' Read HDF5 files
#' @param file path to the HDF5 file to be read
#' @param name name of the attr
#' @param ... parameters passed to \code{\link[rhdf5]{h5read}}
#' @import dplyr
#' @return an array of class \code{brain}
#' @export
#' @examples
#' file <- "~/Data/barresi/AT_1_Probabilities.h5"
#' raw <- read_h5(file)
#' class(raw)
#' dim(raw)
read_h5 <- function(file, name = NULL, ...) {
  if (requireNamespace("hdf5r", quietly = TRUE)) {
    objs <- hdf5r::H5File$new(file, mode = "r")
    if (is.null(name)) {
      # take the first HDF5 dataset
      name <- names(objs)[1]
    }
    x <- objs[[name]]
    if (length(x$dims) > 3) {
      # take the second image
      out <- x[2,,,]
    }

    class(out) <- append("brain", class(out))
    return(out)
  }

  if (requireNamespace("rhdf5", quietly = TRUE)) {
    objs <- rhdf5::h5ls(file) %>%
      dplyr::filter(otype == "H5I_DATASET")
    if (is.null(name)) {
      # take the first HDF5 dataset
      name <- objs$name[1]
    }
    x <- rhdf5::h5read(file, name)
    if (length(dim(x)) > 3) {
      # take the second image
      out <- x[2,,,]
    }

    class(out) <- append("brain", class(out))
    return(out)
  }

  stop("Could not read file. Do you have a HDF5 library available?")
}



#' Tidy 3D brain image data
#' @inheritParams broom::tidy.lm
#' @param threshold probability below which points will be discarded
#' @importFrom tibble as_tibble
#' @importFrom generics tidy
#' @export
#' @examples
#' file <- "~/Data/barresi/AT_1_Probabilities.h5"
#' \dontrun{
#' if (require(dplyr)) {
#'   tidy_brain <- file %>%
#'     read_h5() %>%
#'     tidy()
#'   class(tidy_brain)
#'   glimpse(tidy_brain)
#' }
#' }

tidy.brain <- function(x, type="wildtype", threshold.n=0.9, ...) {
  if (type=="youtoo"){
    if(threshold.n != 0.9){
      threshold = threshold.n
      message(paste0("Set threshold: ", threshold))
    }else{
      threshold = 0.5
      message("Default You-Too threshold: 0.5")
    }
  }else{
    if(threshold.n != 0.9) {  #if threshold is not default, use the input threshold
      threshold = threshold.n
      message(paste0("Set threshold: ", threshold))
    }else{
      threshold = 0.9
      message(("Default threshold: 0.9"))
    }
  }
  res <- x %>%
    as.data.frame.table() %>%
    mutate(x = as.integer(Var1) * 0.13,
           y = as.integer(Var2) * 0.13,
           z = as.integer(Var3) * 0.21) %>%
    select(x, y, z, Freq) %>%
    tibble::as_tibble() %>%
    mutate(gray_val = grDevices::gray(1 - Freq)) %>%
    filter(Freq > threshold)
  class(res) <- append("tbl_brain", class(res))
  return(res)
}

#' Plot a slice of a 3D image
#' @inheritParams graphics::image
#' @param z slice of the 3D image to display
#' @importFrom graphics image
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
  graphics::image(x[,,z])
}

#' Plot a 2D projection of a 3D image
#' @inheritParams plot3d.tbl_brain
#' @param plane a character vector of length 2 indicating the
#' plane to project to (x, y, z)
#' @param show_max Plot the maximum frequency or the average?
#' @export
#' @examples
#' file <- "~/Data/barresi/wildtype/AT/AT_1_Probabilities.h5"
#' \dontrun{
#' tidy_brain <- read_h5(file) %>%
#'   tidy()
#' plot2d(tidy_brain, show_max = TRUE)
#' plot2d(tidy_brain)
#' }

plot2d <- function(x, title="Sample", show_max = FALSE, ...) UseMethod("plot2d")

#' @export
#' @rdname plot2d

plot2d.brain <- function(x, title="Sample", show_max = FALSE, ...) {
  tidy_brain <- tidy(x)
  plot2d(tidy_brain, show_max = FALSE, ...)
}

#' @export
#' @importFrom gridExtra grid.arrange
#' @rdname plot2d
plot2d.tbl_brain <- function(x, title="Sample", show_max = FALSE, ...) {
  gridExtra::grid.arrange(
    plot2d_plane(x, plane = c("x", "z"), show_max, ...),
    plot2d_plane(x, plane = c("x", "y"), show_max, ...),
    plot2d_plane(x, plane = c("y", "z"), show_max, ...),
    nrow = 1,
    top = title
  )
}

#' @export
#' @importFrom rlang !!
#' @import ggplot2
#' @rdname plot2d
plot2d_plane <- function(x, plane = c("x", "z"), show_max = FALSE, ...) {
  depth <- setdiff(c("x", "y", "z"), plane)
  depth_var <- rlang::sym(depth)

  gg_data <- x %>%
    group_by(.dots = plane) %>%
    summarize(N = n(), min_freq = min(Freq),
              avg_depth = mean(!! depth_var, na.rm = TRUE),
              max_depth = max(!! depth_var, na.rm = TRUE))

  if (show_max) {
    plot_var <- "max_depth"
  } else
    plot_var <- "avg_depth"
  labels <- data.frame(var = c("x", "y", "z"),
                       label = c("Lateral (x, microns)",
                                 "Anterior/Posterior (y, microns)",
                                 "Dorsal/Ventral (z, microns)")) %>%
    filter(var %in% plane)
  titles <- data.frame(depth_var = c("x", "y", "z"),
                       title = c("Side View", "Front View", "Top View"))


  plot <- ggplot(gg_data, aes_string(x = plane[1], y = plane[2], color = plot_var), ...) +
    geom_point(alpha = 0.1, size = 0.5) +
    annotate("text", x = 0, y = 0, label = "origin", size = 5) +
    annotate("text", x = 0, y = 0,
             label = paste0("min_prob: ", round(min(pull(ungroup(gg_data), "min_freq")), 2)))+
    scale_color_continuous(guide = FALSE) +
    scale_x_continuous(filter(labels, var == plane[1])$label) +
    scale_y_continuous(filter(labels, var == plane[2])$label) +
    ggtitle(filter(titles, depth_var == depth)$title)
  #+ coord_fixed()

  if (plane == c("x", "y")){
   plot + geom_smooth(method = "lm", formula = y ~ I(x^2) + x, color = "red")
  } else{
   plot + geom_smooth(method = "lm", color = "red")
  }
}

#' Plot a 3D image of a brain
#' @inheritParams rgl::plot3d
#' @export
#' @examples
#' file <- "~/Data/barresi/AT_1_Probabilities.h5"
#' \dontrun{
#' if (require(dplyr)) {
#'   tidy_brain <- file %>%
#'     read_h5() %>%
#'     tidy()
#'   plot3d(tidy_brain)
#' }
#' }

plot3d.tbl_brain <- function(x, ...) {
  x_df <- as.data.frame(x)
  n <- nrow(x_df)
  message(paste("Will now plot", n, "points..."))
  if (n > 250000) {
    warning(paste("You are about to plot", n,
                  "points. Consider increasing the threshold parameter."))
  }
  rgl::plot3d(x_df, type = "p", alpha = x_df$Freq,
              #         width = 900, height = 800,
              xlab = "x", ylab = "y", zlab = "z",
              size = 2, col = x_df$gray_val, ...)
  #  plot3d_model(x_df)
}

#' @rdname plot3d.tbl_brain
#' @export
#' @examples
#' file <- "~/Data/barresi/AT_1_Probabilities.h5"
#' \dontrun{
#' brain <- read_h5(file)
#' plot3d(brain)
#' }


plot3d.brain <- function(x, ...) {
  tidy_brain <- tidy(x)
  return(plot3d(tidy_brain, ...))
}

#' @importFrom mosaic makeFun
#' @importFrom rgl plot3d
#' @importFrom stats lm coef

plot3d_model <- function(xyz, ...) {
  mod_list <- get_jawbone(xyz)
  # quadratic plane
  plot3d(mod_list[[1]], alpha = 0.5, col = "dodgerblue",
         xlim = range(xyz$x), ylim = range(xyz$y), zlim = range(xyz$z),
         add = TRUE)

  # flat plane
  plot3d(mod_list[[2]], alpha = 0.5, col = "green",
         xlim = range(xyz$x), ylim = range(xyz$y), zlim = range(xyz$z),
         add = TRUE)

  # their intersection
  t <- seq(0, max(xyz$x))
  plot3d(mod_list[[3]](t), alpha = 0.5, col = "red", size = 5,
         xlim = range(xyz$x), ylim = range(xyz$y), zlim = range(xyz$z),
         add = TRUE)
  return(mod_list)
}

#' Compute the jawbone model
#' @param xyz a \code{\link{tbl_brain}}
#' @param ... currently ignored
#' @export
#' @examples
#' file <- "~/Data/barresi/AT_1_Probabilities.h5"
#' \dontrun{
#' brain <- read_h5(file)
#' jbone <- get_jawbone(tidy(brain))
#' }

get_jawbone <- function(xyz, ...) {
  # quadratic plane
  mod1 <- stats::lm(z ~ x + y + I(x^2), data = xyz)
  bent_plane <- mosaic::makeFun(mod1)

  # flat plane
  mod2 <- stats::lm(y ~ x + z, data = xyz)
  flat_plane <- mosaic::makeFun(mod2)

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
  f <- function(t) {
    cbind(x = t,
          y = a_prime*t^2 + b_prime*t + c_prime,
          z = (a + c*a_prime)*t^2 + (b + c*b_prime)*t + c*c_prime + d)
  }
  return(list(bent_plane, flat_plane, f))
}

#' @title Reorientation: applying PCA (Principal Component Analysis)
#' @description Principal component analysis (PCA) is a statistical procedure that uses an
#' orthogonal transformation to convert a set of observations of possibly correlated variables
#' (entities each of which takes on various numerical values) into a set of values of linearly
#' uncorrelated variables called principal components.\url{https://en.wikipedia.org/wiki/Principal_component_analysis}

#' @param x a \code{\link{tbl_brain}}
#' @param correctionA identify PCA alignment error type A and make correction \link[cranium]{is_errorA}
#' @param ... arguments passed to \code{\link[stats]{prcomp}}
#' @importFrom tibble as_tibble
#' @export
#' @examples
#' tidy_brain <- read_h5(download_brains(tempdir(), pattern = "101"))%>%
#'        tidy(type="wildtype", threshold = 0.9)
#' ro_tidy_brin <- reorient(tidy_brain)


reorient <- function(x, correctionA = FALSE, ...) {
  pca <- stats::prcomp(~ x + y + z, data = x, scale = FALSE, ...)
  out <- pca$x %>%
    tibble::as_tibble() %>%
    dplyr::rename(x = PC1, y = PC2, z = PC3)
  # fit quadratic model in the xy-plane
  mod <- stats::lm(y ~ x + I(x^2), data = out)
  a <- coef(mod)["I(x^2)"]
  b <- coef(mod)["x"]
  c <- coef(mod)["(Intercept)"]
  # translate to put vertex at origin
  # https://en.wikipedia.org/wiki/Parabola#Parabola_as_graph_of_a_function
  vertex <- c(-b / (2 * a), (4 * a * c - b^2) / (4 * a))
  # translate to center z on commissure
  z_mean <- out %>%
    filter(abs(x) < 50) %>%
    summarize(z_mean = mean(z))
  if (correctionA==FALSE){
    out <- out %>%
      mutate(x = x - vertex[1], y = y - vertex[2], z = z - z_mean$z_mean)
  }else{
    out <- out %>%
      mutate(x = x - vertex[1], y =  vertex[2] - y , z = z - z_mean$z_mean)
  }
  out[, c("Freq", "gray_val")] <- x[, c("Freq", "gray_val")]
  class(out) <- append("tbl_brain", class(out))
  # recompute model after translation
  #xz plane
  attr(out, "quad_mod_xy") <- stats::lm(y ~ x + I(x^2), data = out)
  attr(out, "linear_mod_xy") <- stats::lm(y ~ x, data=out)
  #xz plane
  attr(out, "quad_mod_xz") <- stats::lm(z ~ x + I(x^2), data = out)
  attr(out, "linear_mod_xz") <- stats::lm(z ~ x, data=out)
  #yz plane
  attr(out, "quad_mod_yz") <- stats::lm(z ~ y + I(y^2), data = out)
  attr(out, "linear_mod_yz") <- stats::lm(z ~ y, data=out)
  return(out)
}

#' Change to parabolic coordinates
#' @param obj an object
#' @param ... currently ignored
#' @export
#' @examples
#'
#' file <- "~/Data/barresi/AT_1_Probabilities.h5"
#' \dontrun{
#' if (require(dplyr)) {
#'   tidy_brain <- file %>%
#'     read_h5() %>%
#'     tidy()
#'   plot3d(tidy_brain)
#'   tidy_flat <- reorient(tidy_brain)
#'   plot3d(tidy_flat)
#'   spherical <- change_coordinates(tidy_flat)
#'   summary(spherical)
#'   }
#' }

change_coordinates <- function(obj, ...) {
  a <- coef(attr(obj, "quad_mod"))["I(x^2)"]
  integrand <- function(x, a) { sqrt(1 + (2*a*x)^2); }
  alpha <- lapply(obj$x, stats::integrate, f = integrand, lower = 0, a = a) %>%
    unlist()
  obj$alpha <- alpha[seq(from = 1, to = length(alpha), by = 5)]

}

#' @title Quadratic model evaluation
#' @description Fit quadratic model y=x^2 + x to the data and evaluate how well the model fit the data using adjusted
#' R-Squared, RMSE, number of signals being captured
#' @param data a brain image or a tidy brain
#' @param type wild type or mutant type of zebrafish embryo brain
#' @param threshold.n threshld for signal frequency. Any signal with frequency below threshold value is exluded
#' @export
#' @examples
#' brain <- read_h5(download_brains(tempdir(), pattern = "101"))
#' qmodel(brain)
#'
#' #If you have already tidied your brain
#' tidy_brain <- tidy(brain, type = "wildtype", threshold.n = 0.9)
#' qmodel(tidy_brain)

qmodel <- function(data, type="wildtype", threshold.n=0.9) UseMethod("qmodel")

#' @export
#' @rdname qmodel
qmodel.brain<- function(data, type="wildtype", threshold.n=0.9){
  tidy_brain<-data%>%
    tidy(type, threshold = threshold.n)
  ro_tidy_brain <- data%>%
    tidy(type, threshold = threshold.n)%>%
    reorient()
  model <- lm(y ~ I(x^2) + x, data = tidy_brain) #model applied on original data
  ro_model <- attr(ro_tidy_brain, "quad_mod")  #model applied on reoriented data
  mean.fitted <- mean(ro_model$fitted.values)
  #sigma=sjstats::rse(ro_model). RSE: residual standard error
  sum_tb<- ro_model %>%
    broom::glance() %>%
    mutate(num.signal = nrow(tidy_brain),
           adj.r.squared.original = round(summary(model)$adj.r.squared, 3),#original r square
           quad.coeff = round((summary(ro_model)$coefficients["I(x^2)", "Estimate"]), 4),
           RMSE = sjstats::rmse(ro_model) )
  return(sum_tb)
}

#' @export
#' @importFrom gridExtra grid.arrange
#' @rdname qmodel
qmodel.tbl_brain<- function(tidy_brain, type="wildtype", threshold.n=0.9){
  ro_tidy_brain <- tidy_brain%>%
    reorient()
  model <- lm(y ~ I(x^2) + x, data = tidy_brain) #model applied on original data
  ro_model <- attr(ro_tidy_brain, "quad_mod")  #model applied on reoriented data
  mean.fitted <- mean(ro_model$fitted.values)
  #sigma=sjstats::rse(ro_model). RSE: residual standard error
  sum_tb<- ro_model %>%
    broom::glance() %>%
    mutate(num.signal = nrow(tidy_brain),
           adj.r.squared.original = round(summary(model)$adj.r.squared, 3),#original r square
           quad.coeff = round((summary(ro_model)$coefficients["I(x^2)", "Estimate"]), 4),
           RMSE = sjstats::rmse(ro_model) )
  return(sum_tb)
}


#' @title PCA alignment error type A identification and correction
#' @description Identify Error Type A that the y-axis is flipped y-axis and make correction
#' @param data a brain image
#' @param type wild type or mutant type of zebrafish embryo brain
#' @param threshold.n threshld for signal frequency. Any signal with frequency below threshold value is exluded
#' @export
#' @examples
#' brain <- read_h5(download_brains(tempdir(), pattern = "101"))
#' is_errorA(brain)
#' corrected_brain <- correct_errorA(brain)

is_errorA <- function(data, type="wildtype", threshold.n=0.9) UseMethod("is_errorA")

#' @export
#' @rdname is_errorA
is_errorA.brain<- function(data, type="wildtype", threshold.n=0.9){
  ro_tidy_brain <- data%>%
    tidy(type, threshold = threshold.n)%>%
    reorient()
  ro_model <- attr(ro_tidy_brain, "quad_mod_xy")  #model applied on reoriented data
  quad.coeff = round((summary(ro_model)$coefficients["I(x^2)", "Estimate"]), 4)

  if (quad.coeff < 0) {
    message("Y Axis is flipped")
    return(TRUE)
  } else{
    message("Correct alignment")
    return(FALSE)
  }
}


#' @export
#' @rdname is_errorA
correct_errorA <- function(data, type="wildtype", threshold.n=0.9) UseMethod("correct_errorA")

#' @export
#' @rdname is_errorA
correct_errorA.brain<- function(data, type="wildtype", threshold.n=0.9){
  c_ro_tidy_brain <- data%>%
    tidy(type, threshold = threshold.n)%>%
    reorient(correctionA = TRUE)
  ro_model <- attr(c_ro_tidy_brain, "quad_mod_xy")
  quad.coeff = round((summary(ro_model)$coefficients["I(x^2)", "Estimate"]), 4)
  if (quad.coeff <0){
    message("Y Axis is flipped")
  }else{
    message("Correct alignment")
  }
  return(c_ro_tidy_brain)
}




#' @title Identify Curve in yz plane and xz plane
#' @description Identify Curve in yz plane and xz plane
#' @param data a brain image
#' @param type wild type or mutant type of zebrafish embryo brain
#' @param threshold.n threshld for signal frequency. Any signal with frequency below threshold value is exluded
#' @export
#' @examples

#' @export
#' @rdname is_z_curve
is_z_curve<- function(data, type="wildtype", threshold.n=0.9) UseMethod("is_z_curve")

# @export
# @rdname is_z_curve
#'
#is_z_curve.brain <- function(data, type="wildtype", threshold.n=0.9){
 # ro_tidy_brain <- data%>%
 #   tidy(threshold=threshold.n)
 # sum_tb <- is_z_curve(ro_tidy_brain, type = "wildtype", threshold = threshold.n)
 # return(sum_tb)
#}

#' @export
#' @rdname is_z_curve
is_z_curve.tbl_brain <- function(data, type="wildtype", threshold.n=0.9){
  ro_tidy_brain <- data %>%
    reorient()
  z_mean <- ro_tidy_brain %>%
    skim(z)%>%
    select(stat, formatted)%>%
    spread(stat, formatted)%>%
    pull(mean)

  #Fit models after reorient
  #xz plane
  lm_model_xz <- attr(ro_tidy_brain, "linear_mod_xz") #linear model y = x
  quad_model_xz <- attr(ro_tidy_brain, "quad_mod_xz") #quadratic model y = x^2 + x
  #xy plane
  lm_model_xy <- attr(ro_tidy_brain, "linear_mod_xy")
  quad_model_xy <- attr(ro_tidy_brain, "quad_mod_xy")
  #yz plane
  lm_model_yz <- attr(ro_tidy_brain, "linear_mod_yz")
  quad_model_yz <- attr(ro_tidy_brain, "quad_mod_yz")


  sum_tb <- data.frame(threshold = threshold.n)%>%
    mutate(z_mean = z_mean)%>%
    #xz plane
    mutate(lm.intercept_xz = round((summary(lm_model_xz)$coefficients["(Intercept)", "Estimate"]), 4)) %>%
    mutate(lm.slope_xz = round((summary(lm_model_xz)$coefficients["x", "Estimate"]), 10))%>%
    mutate(quad.slope_xz = round((summary(quad_model_xz)$coefficients["I(x^2)", "Estimate"]), 5)) %>%

    #xy plane
    mutate(lm.intercept_xy= round((summary(lm_model_xy)$coefficients["(Intercept)", "Estimate"]), 4)) %>%
    mutate(lm.slope_xy = round((summary(lm_model_xy)$coefficients["x", "Estimate"]), 10))%>%
    mutate(quad.slope_xy = round((summary(quad_model_xy)$coefficients["I(x^2)", "Estimate"]), 5)) %>%

    #yz plane
    mutate(lm.intercept_yz = round((summary(lm_model_yz)$coefficients["(Intercept)", "Estimate"]), 4)) %>%
    mutate(lm.slope_yz = round((summary(lm_model_yz)$coefficients["y", "Estimate"]), 10))%>%
    mutate(quad.slope_yz = round((summary(quad_model_yz)$coefficients["I(y^2)", "Estimate"]), 5))

  return(sum_tb)
}



#' @title Identify Error B: Curve in yz plane
#' @description Identify Curve in yz plane
#' @param data a brain image
#' @param type wild type or mutant type of zebrafish embryo brain
#' @param threshold.n threshld for signal frequency. Any signal with frequency below threshold value is exluded
#' @export
#' @examples


#' @export
#' @rdname is_errorB
is_errorB <- function(data, type="wildtype", threshold.n=0.9) UseMethod("is_errorB")

#' @export
#' @rdname is_errorB
is_errorB.brain <- function(data, type="wildtype", threshold.n=0.9){
  ro_tidy_brain <- data%>%
    tidy(type, threshold = threshold.n)%>%
    reorient()
  quad_model_xz <- attr(ro_tidy_brain, "quad_mod_xz") #quadratic model y = x^2 + x
  quad.slope_xz = round((summary(quad_model_xz)$coefficients["I(x^2)", "Estimate"]), 5)

  if (quad.slope_xz < -0.00202 | quad.slope_xz > 0.00182) {
    message("Curve in xz plane")
    return(TRUE)
  } else{
    message("Correct alignment")
    return(FALSE)
  }
}

#' @export
#' @rdname is_errorB
is_errorB.tbl_brain <- function(data, type="wildtype", threshold.n=0.9){
  ro_tidy_brain <- data%>%
    reorient()
  quad_model_xz <- attr(ro_tidy_brain, "quad_mod_xz") #quadratic model y = x^2 + x
  quad.slope_xz = round((summary(quad_model_xz)$coefficients["I(x^2)", "Estimate"]), 5)

  if (quad.slope_xz < -0.00202 | quad.slope_xz > 0.00182) {
    message("Curve in xz plane")
    return(TRUE)
  } else{
    message("Correct alignment")
    return(FALSE)
  }
}

