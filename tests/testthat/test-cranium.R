context("test-cranium")

test_that("reading works", {
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
