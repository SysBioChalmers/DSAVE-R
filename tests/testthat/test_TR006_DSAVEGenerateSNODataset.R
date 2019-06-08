test_that("TR006 - DSAVEGenerateSNODataset",{

  ds = LoadAndDownloadBCells();
  ds2 = ds[, 1:100];

  ds3 = DSAVEGenerateSNODataset(ds2);

  UMIs1 = as.vector(colSums(as.matrix(ds2)));
  UMIs2 = as.vector(colSums(as.matrix(ds3)));

  #check that the UMI distribution match
  expect_equal(UMIs1, UMIs2, info = "TR006: UMIDistr mismatch")

  #check that sparse input gives sparse output
  expect_true(is(ds2, 'sparseMatrix') & is(ds3, 'sparseMatrix'), info = "TR006: Failure on sparse input")

  #check that matrix input gives matrix output
  ds4 = as.matrix(ds2);
  ds5 = DSAVEGenerateSNODataset(ds4);
  expect_true(is.matrix(ds4) & is.matrix(ds5), info = "TR006: Failure on non-sparse input")

})
