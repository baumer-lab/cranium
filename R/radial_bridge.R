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

plot2d.brain <- function(x, title="Sample", show_max = FALSE,  ...) {
  tidy_brain <- tidy(x)
  plot2d(tidy_brain, show_max = FALSE, ...)
}

#' @export
#' @importFrom gridExtra grid.arrange
#' @rdname plot2d
plot2d.tbl_brain <- function(x, title="Sample", show_max = FALSE,  ...) {
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
  } else if (show_max == FALSE){
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
      geom_point(aes( x= 0, y = 0), colour="blue", size = 3) +
      #annotate("text", x = 0, y = 0, label = "origin", size = 3) +  #take out text because text is overlapped
      #annotate("text", x = 0, y = min(gg_data$z) * 0.9,
      #         label = paste0("min prob: ", round(min(pull(ungroup(gg_data), "min_freq")), 2)))+
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


#' @title Add attribute of quadratic model
#' @description This function add attribute to tbl_brain data on xy plane, yz plane and xz plane. This facilitate users to
#' access the model information of the data.
#' @param x a \code{\link{tbl_brain}}
#' @export
#' @examples
#' ro_brain_att <- add_attr(ro_brain)
#'
add_attr <- function(x){
  attr(x, "quad_mod_xy") <- stats::lm(y ~ x + I(x^2), data = x)
  attr(x, "linear_mod_xy") <- stats::lm(y ~ x, data= x)
  #xz plane
  attr(x, "quad_mod_xz") <- stats::lm(z ~ x + I(x^2), data =x)
  attr(x, "linear_mod_xz") <- stats::lm(z ~ x, data= x)
  #yz plane
  attr(x, "quad_mod_yz") <- stats::lm(z ~ y + I(y^2), data = x)
  attr(x, "linear_mod_yz") <- stats::lm(z ~ y, data= x)
  return(x)
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


reorient <- function(x,...) {
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
  out <- out %>%
    mutate(x = x - vertex[1], y = y - vertex[2], z = z - z_mean$z_mean)

  out[, c("Freq", "gray_val")] <- x[, c("Freq", "gray_val")]
  class(out) <- append("tbl_brain", class(out))

  # recompute model after translation
  out <- add_attr(out)
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
           quad.coef = round((summary(ro_model)$coefficients["I(x^2)", "Estimate"]), 4),
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
           quad.coef = round((summary(ro_model)$coefficients["I(x^2)", "Estimate"]), 4),
           RMSE = sjstats::rmse(ro_model) )
  return(sum_tb)
}


#' @title PCA alignment error type A identification and correction
#' @description Identify Error Type A that the y-axis is flipped y-axis and make correction. If the object is brain type (raw data), will
#' perform tidy and PCA reorient functions. If the object is a tbl_brain class (already performed PCA reorientation), will not
#' perform PCA reorientation again.
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
is_errorA.brain<- function(data, type="wildtype", threshold.n = 0.9){
  ro_tidy_brain <- data%>%
    tidy(type, threshold = threshold.n)%>%
    reorient()
  ro_model <- attr(ro_tidy_brain, "quad_mod_xy")  #model applied on reoriented data
  quad.coef = round((summary(ro_model)$coefficients["I(x^2)", "Estimate"]), 4)

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
is_errorA.tbl_brain<- function(data, type="wildtype", threshold.n=0.9){
  # fit quadratic model in the xy-plane, yz-plane, and xz-plane
  data <- add_attr(data)

  ro_model <- attr(data, "quad_mod_xy")  #model applied on reoriented data
  quad.coef = round((summary(ro_model)$coefficients["I(x^2)", "Estimate"]), 4)

  if (quad.coef < 0) {
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
  quad.coef = round((summary(ro_model)$coefficients["I(x^2)", "Estimate"]), 4)
  if (quad.coef < 0){
    message("Y Axis is flipped")
  }else{
    message("Correct alignment")
  }
  return(c_ro_tidy_brain)
}


#' @export
#' @rdname is_errorA
correct_errorA <- function(data, type = "wildtype", threshold.n=0.9) UseMethod("correct_errorA")

#' @export
#' @rdname is_errorA
correct_errorA.brain <- function(data, type = "wildtype", threshold.n=0.9){
  ro_tidy_brain <- data%>%
    tidy(type, threshold.n)%>%
    reorient()

  r.data <- correct_errorA.tbl_brain(ro_tidy_brain)
  return(r.data)
}


#' @export
#' @rdname is_errorA
correct_errorA.tbl_brain <- function(data, type = "wildtype", threshold.n=0.9){
  rotation_angle_radian <- pi
  #radian to degree
  rotation_angle_degree <- 180

  matrix.x <- diag(3)
  matrix.x[2,2] <- cos(rotation_angle_radian)
  matrix.x[2,3] <- -sin(rotation_angle_radian)
  matrix.x[3,2] <- sin(rotation_angle_radian)
  matrix.x[3,3] <- cos(rotation_angle_radian)

  #prepare data for matrix multiplication
  data.matrix<-data%>%
    select(x,y,z)%>%
    as.matrix()

  #keep information other than coordinates in seperate dataset
  data_freq_gray<- data%>%
    select(Freq, gray_val)

  r.data <- data.matrix %*% matrix.x %>%
    tibble::as_tibble() %>%
    bind_cols(data_freq_gray)

  names(r.data) <- names(data)

  class(r.data) <- append("tbl_brain", class(r.data))

  return(r.data)
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


#' @title create rotation matrix
#' @description creat 3 by 3 rotation matrix for error B correction
#' @param rotation_angle_radian radian of the estimated rotation angle
#' @param axis the axis that you are rotating around
#' @export
#' @examples
#' degree = 45
#' radian = deg2rad(degree)
#' ro_matrix <- create_ro_matrix(radian, axis = "y")
create_ro_matrix <- function(rotation_angle_radian, axis = "x"){
  if(axis == "x"){
    matrix.x <- diag(3)
    matrix.x[2,2] <- cos(rotation_angle_radian)
    matrix.x[2,3] <- -sin(rotation_angle_radian)
    matrix.x[3,2] <- sin(rotation_angle_radian)
    matrix.x[3,3] <- cos(rotation_angle_radian)
    message("rotate around x axis")
    return(matrix.x)
  } else if (axis == "y"){
    matrix.x <- diag(3)
    matrix.x[1,1] <- cos(rotation_angle_radian)
    matrix.x[1,3] <- sin(rotation_angle_radian)
    matrix.x[3,1] <- -sin(rotation_angle_radian)
    matrix.x[3,3] <- cos(rotation_angle_radian)
    message("rotate around y axis")
    return(matrix.x)
  }else if (axis == "z"){
    matrix.x <- diag(3)
    matrix.x[1,1] <- cos(rotation_angle_radian)
    matrix.x[2,1] <- sin(rotation_angle_radian)
    matrix.x[1,2] <- -sin(rotation_angle_radian)
    matrix.x[2,2] <- cos(rotation_angle_radian)
    message("rotate around z axis")
    return(matrix.x)
  }
}


#' @title Identify Error B: Curve in yz plane
#' @description Identify Curve in yz plane
#' @param data a brain image
#' @param type wild type or mutant type of zebrafish embryo brain
#' @param threshold.n threshld for signal frequency. Any signal with frequency below threshold value is exluded
#' @param range_ref error level.Options are 1. 95% Confidence Interval; 2.Inter-quatile range (IQR); 3. outlier range. Default to
#' outlier, which could be interpreted as: if the quadratic coefficient of the youtoo type sample is within the normal range of
#' wildtype's quadratic coeffcients, we say there is no error. Otherwise, if the quadratic coefficient is considered as an outlier
#' in the wildtype quadratic coefficients, we say there is a an error B in the youtoo type sample.
#' @export
#' @examples
#' is_errorB(data)

#' @export
#' @rdname is_errorB
is_errorB <- function(ro_data, type="wildtype", threshold.n=0.9) UseMethod("is_errorB")

#' @export
#' @rdname is_errorB
is_errorB.brain <- function(ro_data, type="wildtype", threshold.n=0.9){
  ro_tidy_brain <- data%>%
    tidy(type, threshold = threshold.n)%>%
    reorient()
  is_errorB(ro_tidy_brain)
}


#' @export
#' @rdname is_errorB
is_errorB.tbl_brain <- function(ro_data, type="wildtype", threshold.n=0.9){
  ro_data <- add_attr(ro_data)
  quad_model <- attr(ro_data, "quad_mod_xz") #quadratic model y = x^2 + x
  quad.coef= round((summary(quad_model)$coefficients["I(x^2)", "Estimate"]), 5)

  ro_model <- attr(ro_data, "quad_mod_xz")
  quad.coef = round((summary(ro_model)$coefficients["I(x^2)", "Estimate"]), 4)  #a


  if (abs(quad.coef) <= 0.0001){
    message("Alignment is correct. There is no curvature in the xz-plane.")
    return(FALSE)
  }else{
    message("There is a curvature in xz plane")
    return(TRUE)
  }
}


#' @export
#' @rdname is_errorB
correct_errorB.tbl_brain <- function(ro_data, r_angle_mpt=1) UseMethod("correct_errorB")

#' @export
#' @rdname is_errorB
correct_errorB.tbl_brain <- function(ro_data, r_angle_mpt=1){
  ro_data <- add_attr(ro_data)

  #extract quadratic model coefficients from xy plane
  ro_model <- attr(ro_data, "quad_mod_xz")
  quad.coef = round((summary(ro_model)$coefficients["I(x^2)", "Estimate"]), 4)  #a

  if (abs(quad.coef) <= 0.0001){
    message("Before rotation, the quadratic coefficient is 0, no need for further correction. Return original data.")
    return (ro_data)
  }else{
  message("Before rotation, the quadratic coefficient is", " ", quad.coef, ".")

  x.coeff = round((summary(ro_model)$coefficients["x", "Estimate"]), 100)#b
  y_intercept = round((summary(ro_model)$coefficients["(Intercept)", "Estimate"]), 4) #c

  discriminant <- x.coeff^2 - 4* quad.coef * y_intercept
  if(discriminant < 0){
    message("Discriminant <0. Can't make more corrections. Return original data.")
    return(ro_data)
  }else{
    # x intercept
    x_intercept <- (-x.coeff + sqrt(x.coeff^2 - 4* quad.coef * y_intercept))/ (2*quad.coef)

    #slope of the line connecting m and n
    slope <- abs(y_intercept/x_intercept)

    #slope to angle (radian)
    rotation_angle_radian <- r_angle_mpt * atan(slope)  #r_angle_mpt is a tuning parameter, >1 leads to larger rotation angel
    #radian to degree
    #rotation_angle_degree <- rad2deg(rotation_angle_radian)


    if(quad.coef > 0 ){
      message("quadratic concave up, will perform error B correction")
      #rotation matrix
      matrix.x <- create_ro_matrix(rotation_angle_radian, axis = "x")

      data.matrix<- ro_data %>%
        select(x,y,z)%>%
        as.matrix()
      data_freq_gray<- ro_data%>%
        select(Freq,gray_val)
      r.data <- data.matrix %*% matrix.x %>%
        tibble::as_tibble() %>%
        bind_cols(data_freq_gray)
      names(r.data) <- names(ro_data)
      class(r.data) <- append("tbl_brain", class(r.data))
      r.data <- add_attr(r.data)

      #check the quadratic coefficient after rotation
      r.ro_model <- attr(r.data, "quad_mod_xz")
      r.quad.coef = round((summary(r.ro_model)$coefficients["I(x^2)", "Estimate"]), 4)  #a
      message("After rotation, the quadratic coefficient is", " ", r.quad.coef, ".")

      return(r.data)
    }
    else if (quad.coef < 0 ){
      message("quadratic concave down, will need to perform translation before performing error B correction")
      ro_data <- ro_data%>%
        mutate(z = z - y_intercept)  # if y_intercept <0, we translate up, if y intercept >0, we translate down
      rotation_angle_radian <- -rotation_angle_radian
      #rotation matrix
      matrix.x <- create_ro_matrix(rotation_angle_radian, axis = "x")

      data.matrix<- ro_data %>%
        select(x,y,z)%>%
        as.matrix()

      data_freq_gray<- ro_data%>%
        select(Freq,gray_val)

      r.data <- data.matrix %*% matrix.x %>%
        tibble::as_tibble() %>%
        bind_cols(data_freq_gray)
      names(r.data) <- names(ro_data)
      r.data <- r.data%>%
        mutate(z = z + y_intercept)
      class(r.data) <- append("tbl_brain", class(r.data))
      r.data <- add_attr(r.data)
      #check the quadratic coefficient after rotation
      r.ro_model <- attr(r.data, "quad_mod_xz")
      r.quad.coef = round((summary(r.ro_model)$coefficients["I(x^2)", "Estimate"]), 4)  #a
      message("After rotation, the quadratic coefficient is", " ", r.quad.coef, ".")

      return(r.data)
  }
  }
  }
}





#' @title Read in and tidy all raw samples
#' @description Read in and tidy all samples. You can choose to perform PCA reorientation or not. Note that this function is only for raw data (untidied data).
#' Function will output a list containing all samples in the type you decide (wildtype or youtoo), in
#' the format of dataframes. Processed samples "ro_brain_ls_wt.rds" and "ro_brain_ls_yt.rds" are generated from the
#' raw data using this function.
#' @param download TRUE or FALSE, If TRUE, you will download data to your current directory. If FALSE, you will need to specify the directory
#' where you keep you files locally by setting "local_directory".
#' @param local_directory If you have samples in your locally, enter the directory. Note than this should be a folder that keeps the samples
#' of the same sample type such as wild type. You should not keep both wild type and youtoo type in the same folder.
#' @param type wild type or youtoo type of zebrafish embryo brain
#' @param threshold.n threshld for signal frequency. Any signal with frequency below threshold value is exluded.
#' Default to 0.5
#' @param PCA_reorient TRUE or FALSE. If TRUE, we will tidy and then perform PCA reorientation to your data. If FALSE, we will only tidy your data.
#' @export
#' @examples
#' ro_brain_ls_wt <- process_all_sample(download = FALSE, local_directory = NA, type = "wildtype", threshold.n = 0.5, PCA_reorient = TRUE)
#' ro_brain_ls_yt <- process_all_sample(download = FALSE, local_directory = NA, type = "youtoo", threshold.n = 0.5, PCA_reorient = TRUE)
#' saveRDS(ro_brain_ls_wt , file = "ro_brain_ls_wt.rds")
#' saveRDS(ro_brain_ls_yt , file = "ro_brain_ls_yt.rds")
process_all_sample <- function(download = FALSE, local_directory = NA, type = "youtoo", threshold.n = 0.5, PCA_reorient = TRUE){
  #read in files
    if (download == TRUE){
      local_directory <- paste("./sample_download_", type, sep= "")
      dir.create(local_directory)
     # setwd(local_directory)

      if (type == "wildtype"){
        wildtype_files <- download_brains(path = local_directory, pattern = "wildtype")
      }else if (type == "youtoo") {
        message("We have not hosted any youtoo sample yet. Please try download wildtype samples.")
      }
    }else {
      message("We will read in files from your local directory")
    }
    message(paste("Read in and tidy data--", type, sep = ""))
    file_names <- list.files(local_directory, pattern=".h5", full.names=TRUE)
    tidy_brain_ls<- list()
    for (i in (1:length(file_names))) {
      tidy_brain_ls[[i]] <- file_names[i] %>%
        read_h5() %>%
        tidy.brain(threshold.n = threshold.n)     ##note that we should be able to call the global function tidy
    }
    if (PCA_reorient == FALSE){
      return(tidy_brain_ls)  #if do not perform PCA reorientation, return tidy data
    }else{
      print(paste("Process PCA reorientation--", type, sep = ""))
      #PCA reorientation
      ro_brain_ls <- tidy_brain_ls %>%
        lapply(reorient)
      return(ro_brain_ls)
  }
}



#' @title Plot all samples
#' @description Plot all samples together and output as pdf
#' @param x the list of samples you want to plot
#' @param type wild type or youtoo type of zebrafish embryo brain
#' @param file_directory samples' local directory, we will need this to extract the name of each of your sample. This should correspond to the list
#' of sampels.
#' @param save_directory Set a directory for the plots, otherwise, the plot will be saved at the current directory.
#' @export
#' @examples
#' plot_all_sample(tidy_brain_ls, type = "wildtype", file_directory = NA, save_directory = NA)
#' plot_all_sample(tidy_brain_ls, type = "youtoo", file_directory = NA, save_directory = NA)
#'
plot_all_sample <- function(x = tidy_brain_ls, type = "wildtype", file_directory = NA, save_directory = NA){
  #open sample_plot folder in cranium
  if (is.na(file_directory)){
    message("Data unfound. If you not have sample data locally, please download the data first.")
  }else {
    message(paste("Data found in directory:", file_directory))
  }

  #create another folder in sample_plot folder with timestamp
  #everytime runing this function will create a folder to store the plots
  time <- Sys.time()
  if (is.na(save_directory)){
    new_directory <- paste("./plot_",type, "_", time, sep="")
  }else{
    new_directory <- paste(save_directory, "/plot_",type, "_", time, sep="")
  }
  dir.create(new_directory)

  file_names <- list.files(file_directory, pattern=".h5", full.names=FALSE)

  pdf(paste(new_directory, "/",type, ".pdf", sep=""), height = 4, width = 15.5)
    for (i in (1:length(file_names))){
      plot2d(x[[i]], title=paste("Sample ",type, " ", file_names[i], sep=""))
    }
  dev.off()
}


#' @title Pipeline for error A, error B correction
#' @description Detect error A and error B sequentially and make corresponding corrections. Return corected brain data.
#' @param x brain data after tidied and PCA reoriented.
#' @param B.correct.n number of error B correction being applied. Since error B correction requires several iterations,
#' suggested number of iteration is 8. Defalt to 8 iterations.
#' @export
#' @examples
#'
#'
pipeline_correctionAB <- function(x = ro_brain_yt, B.correct.n = 8){
  #check for error A
  if (is_errorA(x)){
    message("Detect error A")
    #make correction for error A
    x_step1 <- correct_errorA(x, threshold.n=0.5)
    #confirm correction is made
    if (is_errorA(x_step1)){
      message("Error A correction failed. Return brain data before processing.")
      return(x)
    }else {
      message("Error A correction succeed. Data is corrected for error A.")
      x_step1 <- x_step1
    }
  }else{
    #no error A, keep the original brain data
    message("Error A not detected")
    x_step1 <- x
  }

  #check for error B
  if (is_errorB(x_step1)){
    # there is errorB
    # make correction for error B
    x.list <- list()
    x.list[[1]] <- x_step1
    for (i in 2 : B.correct.n){
      x.list[[i]] <-  correct_errorB.tbl_brain(x.list[[i-1]])
    }
    x_step2 <- x.list[[B.correct.n]]
    if (is_errorB(x_step2)){
      message(paste("After performing errorB correction", B.correct.n, "times, error B is still detected, please increase the the B.correct.n parameter. Return brain data in the previous step."))
      return(x_step2)
    }else{
      message(paste("After performing errorB correction", B.correct.n, "times, we successfully corrected brain data for error B. Return corrected brain data."))
      return(x_step2)
    }
  }else{
    # there is no errorB
    message("Error B is not detected, return original brain data.")
    return(x_step1)
  }
}
