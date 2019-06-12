test_that("TR008 - DSAVEGetSingleCellDivergence",{

  d <- rbind(c(0,1,2,1),c(3,3,1,2));
  ds <- as.matrix(d);

  lls = DSAVEGetSingleCellDivergence(ds, 3, silent = TRUE);

  ds11 <- matrix(1,30000,2);
  lls11 = DSAVEGetSingleCellDivergence(ds11, 400, silent = TRUE);


  d2 <- rbind(c(0,14,2,1),c(17,0,1,2));
  ds2 <- as.matrix(d2);

  lls2 = DSAVEGetSingleCellDivergence(ds2, 3, silent = TRUE);

  d3 <- rbind(c(0,3,2,1),c(3,0,1,2));
  ds3 <- as.matrix(d3);

  lls3 = DSAVEGetSingleCellDivergence(ds3, 3, silent = TRUE);

  d4 <- rbind(c(0,14,5,1),c(17,0,1,2));
  ds4 <- as.matrix(d4);
  lls4 = DSAVEGetSingleCellDivergence(ds4, 6, silent = TRUE);

  #check that more divergent cells get higher divergence
  expect_true((lls[4] > lls[3]) & (lls[4] > lls[1]) &&
              (lls[4] >= lls[2]), info = "TR008: Value not ok")

  #Also explicitly check the values against matlab.
  expVals = c(-1.1240803,-0.8139254,-1.6023828,-0.8139254);
  expect_equal(lls, expVals, info = "TR008: Values not same as in matlab", tolerance=1e-4)

  #check that downsampling gives equal results for the case set up
  expect_equal(lls2, lls3, info = "TR008: Down-sampling not ok", tolerance=1e-4)

  #check that NA is given for only the cells with too few counts
  expect_true(is.na(lls4[4]) & !is.na(lls4[3]), info = "TR008: NA not ok")


  #also make sure it runs for a real dataset
  ds10 = loadOrDownloadB10k();
  ds10 = ds10[, 1:100];
  lls10 = DSAVEGetSingleCellDivergence(ds10, 200, silent = TRUE);
  #check that we have no NA
  expect_equal(sum(is.na(lls10)), 0, info = "TR008: NA Values in b10k")

  })
