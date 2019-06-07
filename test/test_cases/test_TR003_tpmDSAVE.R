# Create a very simple dataset with only two cells and one valid gene, this should be
# deterministic with the same value for all points! This does not fully
# test the function, but at least partly.
test_that("TR003 - TPM",{
  d <- rbind(c(1,2,3,4),c(5,6,7,8), c(9,10,0,12), c(5,32,0,76));
  s1 <- as.matrix(d);
  s2 <- tpmDSAVE(s1);

  expect_equal(s2[2,1],250000, info = "TR003: TPM", tolerance=1e-20)
  expect_equal(s2[3,4],120000, info = "TR003: TPM", tolerance=1e-20)

})
