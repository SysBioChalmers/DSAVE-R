test_that("TR004 - DSAVEGetGeneVariation",{
  #Create a dataset where we can easily calculate the p value
  #and see that it gets somewhat right.
  #There is a risk that this test will have problems with stability, since
  #there are random numbers involved

  d <- rbind(c(0,4,0,0,0),c(0,0,0,0,0),c(4,0,4,4,4));
  ds <- as.matrix(d);
  row.names(ds) = c('A','B','C')


  #the number of combinations in total are 5^4 = 625
  #the number of combinations with 4 in the same pile is 5
  #the p value should thus be 5/625 = 0.0080
  results = DSAVEGetGeneVariation(ds,0,100000,10000);
  expGenes = c('A','C')
  #calculate the expected CV
  CVDS = log(sd(ds[1,] * 250000) / (200000 + 0.05) + 1);
  #calculate the expected SNO CV mean:

  #generate all combinations (not all are unique!) using a loop; if all combinations are there once
  #it will be the true mean
  avgSNO = matrix(0,625,5)
  #the loop variables represents the cell index where each of the four counts should be
  #placed
  for (i in 1:5) {
    for (j in 1:5) {
      for (k in 1:5) {
        for (m in 1:5) {
          index = (i-1)*125 + (j-1)*25 + (k-1)*5 + m;
          avgSNO[index,i] = avgSNO[index,i] + 1;
          avgSNO[index,j] = avgSNO[index,j] + 1;
          avgSNO[index,k] = avgSNO[index,k] + 1;
          avgSNO[index,m] = avgSNO[index,m] + 1;
        }
      }
    }
  }

  CVSNO = log(apply(avgSNO * 250000, 1, sd) / (200000 + 0.05) + 1);
  CVDiffExp = CVDS - mean(CVSNO);
  #now test with different number of UMIs
  #here, the result should be 1/1000 = 0.001
  d2 <- rbind(c(0,3,0,0,0),c(0,0,0,0,0),c(6,0,6,6,9));
  ds2 <- as.matrix(d2);
  row.names(ds2) = c('A','B','C')


  results2 = DSAVEGetGeneVariation(ds2,0,100000,10000);

  expect_equal(results$pVals[1], 0.008, info = "TR004: pVals same UMI count", tolerance=2e-3)
  expect_equal(results2$pVals[1], 0.001, info = "TR004: pVals different UMI count", tolerance=0.0007)
  expect_equal(results$genes, expGenes, info = "TR004: gene mismatch")
  expect_equal(as.numeric(results$logCVDifference[1]), CVDiffExp, info = "TR004: logCV", tolerance=0.01)

})
