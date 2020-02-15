test_that('datatypes', {
    # processFCS
    expect_type(
      processFCS(files, assayname, colsDiscard,
        colsRetain, newColnames),
      c('character'))
    expect_type(processFCS(filter, transformation),
      c('logical'))
    expect_type(processFCS(transFun),
      c('closure'))
    expect_type(
      processFCS(bgNoiseThreshold, euclideanNormThreshold,
        asinhFactor, downsample),
      c('integer'))
    expect_type(processFCS(metadata),
      c('S4', 'list'))
    expect_type(processFCS(downsampleVar),
      c('double'))
    expect_gt(
      processFCS(bgNoiseThreshold, euclideanNormThreshold,
        asinhFactor, downsample, downsampleVar),
      0)

    # importData
    expect_type(processFCS(mat, metadata),
      c('S4', 'list'))
    expect_type(processFCS(mat, downsampleVar),
      c('double'))
    expect_type(processFCS(assayname),
      c('character'))
    expect_gt(processFCS(downsampleVar), 0)
})
