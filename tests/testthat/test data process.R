require(testthat)
require(MSstatsQC)

test_that("camelCaseSplit splits the camel Case correctly",{
  expect_equal(camelCaseSplit("HelloWorld"),"Hello world")
}
)

test_that("punc_remove removes all punctuations from a string",{
  expect_equal(punc_remove("Best_RT"),"Best RT")
}
)

test_that("clearString changes all upper level letters to lower
          case and removes the punctuations",{
            expect_equal(clearString("thisIs_MSstatsQC.Package"),
                         "this is msstats qc package")
}
)

test_that("guessColumnName is returning the column names that we want",{
  expect_equal(guessColumnName("prucurs"),"Precursor")
  expect_equal(guessColumnName("minsttime"),"MinStartTime")
  expect_equal(guessColumnName("aqired.time"),"AcquiredTime")
  expect_equal(guessColumnName("max.time"),"MaxEndTime")
  expect_equal(guessColumnName("wdrft"),"wdrft")
}
)
