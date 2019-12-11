# Create a very simple dataset with only two cells and one valid gene, this should be
# deterministic with the same value for all points! This does not fully
# test the function, but at least partly.
test_that("TR002 - Total variation in single-cell data",{
  d <- rbind(c(1,2),c(999998.9,999997.9), c(0.1,0.1));
  ds <- as.matrix(d);
  rownames(ds) = c('A','B','C');


  val = DSAVEGetTotalVariationPoolSize(ds, 100, 0.5, 1, silent = T)$Rs;
  expVal = log(2.05/1.05);
  # compensate for round off effects:
  #expVals = round(expVals,10);
  #vals = round(vals,10);
  expect_equal(val,expVal, info = "TR002: DSAVEGetTotalVariationPoolSize", tolerance=1e-10)

})
