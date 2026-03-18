test_that("Peak_RankR validates tsv_file_df type", {
  expect_error(
    Peak_RankR(tsv_file_df = "not_a_df", bw_table = data.frame()),
    "`tsv_file_df` must be a data frame"
  )
})

test_that("Peak_RankR catches missing required columns", {
  bad_df <- data.frame(chr = "chr1", start = 100, end = 200)
  bw     <- data.frame(file_path = "a.bw", sample_id = "A")
  expect_error(
    Peak_RankR(tsv_file_df = bad_df, bw_table = bw,
                group_by_column_name = "cell_type"),
    "missing from `tsv_file_df`"
  )
})

test_that("Peak_RankR catches missing bw_table columns", {
  peaks  <- data.frame(chr = "chr1", start = 100, end = 200,
                       cell_type = "A")
  bad_bw <- data.frame(path = "a.bw", id = "A")  # wrong names
  expect_error(
    Peak_RankR(tsv_file_df = peaks, bw_table = bad_bw,
                group_by_column_name = "cell_type"),
    "missing from `bw_table`"
  )
})

test_that("Peak_RankR catches invalid weights — wrong length", {
  peaks <- data.frame(chr = "chr1", start = 100, end = 200,
                      cell_type = "A")
  bw    <- data.frame(file_path = "a.bw", sample_id = "A")
  expect_error(
    Peak_RankR(tsv_file_df = peaks, bw_table = bw,
                group_by_column_name = "cell_type",
                weights = c(1, 1)),
    "length 3"
  )
})

test_that("Peak_RankR catches negative weights", {
  peaks <- data.frame(chr = "chr1", start = 100, end = 200,
                      cell_type = "A")
  bw    <- data.frame(file_path = "a.bw", sample_id = "A")
  expect_error(
    Peak_RankR(tsv_file_df = peaks, bw_table = bw,
                group_by_column_name = "cell_type",
                weights = c(1, -1, 1)),
    "non-negative"
  )
})

test_that("Peak_RankR catches empty tsv_file_df", {
  empty_df <- data.frame(chr = character(), start = integer(),
                         end = integer(),
                         cell_type = character())
  bw <- data.frame(file_path = "a.bw", sample_id = "A")
  expect_error(
    Peak_RankR(tsv_file_df = empty_df, bw_table = bw,
                group_by_column_name = "cell_type"),
    "has no rows"
  )
})

test_that("Peak_RankR accepts custom column names (validation only)", {
  peaks <- data.frame(
    chrom      = "chr1",
    chromStart = 100L,
    chromEnd   = 200L,
    score      = 3.5,
    subclass   = "TypeA"
  )
  bw <- data.frame(file_path = "dummy.bw", sample_id = "TypeA")

  expect_no_error(
    PeakRankR:::.validate_inputs(
      tsv_file_df          = peaks,
      group_by_column_name = "subclass",
      bw_table             = bw,
      weights              = c(1, 1, 1),
      chr_col              = "chrom",
      start_col            = "chromStart",
      end_col              = "chromEnd",
      magnitude_col        = "score"
    )
  )
})

# ── .normalise ────────────────────────────────────────────────────────────────

test_that(".normalise returns values in [0, 1]", {
  x      <- c(1, 5, 3, 10, 2)
  result <- PeakRankR:::.normalise(x)
  expect_true(all(result >= 0 & result <= 1))
  expect_equal(min(result), 0)
  expect_equal(max(result), 1)
})

test_that(".normalise handles constant vector", {
  expect_equal(PeakRankR:::.normalise(c(5, 5, 5)), c(0, 0, 0))
})

test_that(".normalise handles NA values", {
  x      <- c(1, NA, 5)
  result <- PeakRankR:::.normalise(x)
  expect_true(is.na(result[2]))
  expect_equal(result[1], 0)
  expect_equal(result[3], 1)
})

test_that(".normalise handles single-element vector", {
  expect_equal(PeakRankR:::.normalise(42), 0)
})

# ── check_bedtools ────────────────────────────────────────────────────────────

test_that("check_bedtools errors clearly when bedtools is absent", {
  old_path <- Sys.getenv("PATH")
  old_opt  <- getOption("bedtools.path")
  on.exit({
    Sys.setenv(PATH = old_path)
    options(bedtools.path = old_opt)
  })

  Sys.setenv(PATH = "")
  options(bedtools.path = NULL)

  expect_error(check_bedtools(), "bedtools was not found")
  expect_error(check_bedtools(), "https://bedtools.readthedocs.io")
})

test_that("check_bedtools warns for bad options path before falling back", {
  old_opt <- getOption("bedtools.path")
  on.exit(options(bedtools.path = old_opt))
  options(bedtools.path = "/nonexistent/path/")

  expect_warning(
    tryCatch(check_bedtools(), error = function(e) NULL),
    "bedtools not found at the path set in options"
  )
})

# ── Output structure ──────────────────────────────────────────────────────────

test_that("output has expected new columns when bw files are missing", {
  # When no matching bw files exist for a group, Peak_RankR should still
  # return all rows with NA specificity/sensitivity and a magnitude-only rank.
  peaks <- data.frame(
    chr       = c("chr1", "chr1"),
    start     = c(100L, 500L),
    end       = c(200L, 600L),
    cell_type = c("A", "A")
  )
  bw <- data.frame(file_path = "no_match.bw", sample_id = "B")  # no match for A

  # Patch check_bedtools in PeakRankR namespace to avoid system dependency
  old_fn <- get("check_bedtools", envir = asNamespace("PeakRankR"))
  on.exit(
    assignInNamespace("check_bedtools", old_fn, ns = "PeakRankR"),
    add = TRUE
  )
  assignInNamespace(
    "check_bedtools",
    function(...) invisible("/mock/bedtools"),
    ns = "PeakRankR"
  )

  result <- suppressWarnings(
    Peak_RankR(
      tsv_file_df          = peaks,
      group_by_column_name = "cell_type",
      bw_table             = bw,
      verbose              = FALSE
    )
  )

  # Output should have all 5 new columns
  expect_true(all(c("specificity_score", "sensitivity_score",
                    "magnitude_score", "PeakRankR_score",
                    "PeakRankR_rank") %in% names(result)))
  # Row count must be preserved
  expect_equal(nrow(result), nrow(peaks))
  # No bw match → specificity and sensitivity are NA
  expect_true(all(is.na(result$specificity_score)))
  expect_true(all(is.na(result$sensitivity_score)))
  # With no bw match, both magnitude scores are 0 — ranks tied
  expect_true(all(result$PeakRankR_rank >= 1L))
})


# ── Column conflict validation ────────────────────────────────────────────────

test_that("Peak_RankR catches duplicate column arguments", {
  peaks <- data.frame(chr = "chr1", start = 100L, end = 200L,
                      cell_type = "A")
  bw <- data.frame(file_path = "a.bw", sample_id = "A")
  # chr_col and group_by_column_name both set to "chr" — should error
  expect_error(
    Peak_RankR(tsv_file_df = peaks, bw_table = bw,
                group_by_column_name = "chr", chr_col = "chr"),
    "Column name conflict"
  )
})

# ── background_group mismatch warning ─────────────────────────────────────────

test_that("Peak_RankR warns when background_group has unmatched values", {
  peaks <- data.frame(chr = "chr1", start = 100L, end = 200L,
                      cell_type = "A")
  bw <- data.frame(file_path = "no_match.bw", sample_id = "B")

  old_fn <- get("check_bedtools", envir = asNamespace("PeakRankR"))
  on.exit(assignInNamespace("check_bedtools", old_fn, ns = "PeakRankR"),
          add = TRUE)
  assignInNamespace("check_bedtools",
                    function(...) invisible("/mock/bedtools"),
                    ns = "PeakRankR")

  expect_warning(
    suppressWarnings(
      Peak_RankR(tsv_file_df = peaks, bw_table = bw,
                  group_by_column_name = "cell_type",
                  background_group = c("A", "TYPO_GROUP"),
                  verbose = FALSE)
    ),
    "background_group contains values not found"
  )
})
