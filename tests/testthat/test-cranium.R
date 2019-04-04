context("test-cranium")

test_that("reading works", {
  # file <- "https://s3.us-east-2.amazonaws.com/deltascope/AT_04_Probabilities.h5"
  # download.file(file, destfile = tempfile())

  file <- "~/Data/barresi/AT_1_Probabilities.h5"
  if (file.exists(file)) {
    brain <- read_h5(file)
    expect_is(brain, "brain")
    tidy_brain <- tidy(brain)
    expect_is(tidy_brain, "tbl_brain", "tbl_df")
    expect_equal(ncol(tidy_brain), 5)

    expect_is(plot2d(tidy_brain), "gtable")
  }
})
