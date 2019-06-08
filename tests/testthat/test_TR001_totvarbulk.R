test_that("TR001 - Total variation in bulk",{
  d <- rbind(c(1,1,1,1,2,2,2,2),c(300,300,300,300,300,300,300,300), c(0,0,0,0,0,0,0,0));
  s <- as.matrix(d);
  rownames(s) = c('A','B','C');

  t1vs1 = DSAVEGetTotalVariationFromBulk(s, FALSE, 250, 0.5, rescale = FALSE);
  #calculate (each value will be compared with 3 samples of the same value and 4 samples with the double value)
  exp1vs1 = log(2.05/1.05)*4/7;

  #4 vs 4
  t4vs4 = DSAVEGetTotalVariationFromBulk(s, TRUE, 250, 0.5, rescale = FALSE);
  #calculate this differently than in the function, with loops :)
  #calculate the distribution over number of ones (vs twos) in the first of
  #the two sets for each combination
  numones = rep(0,5);#represents 0 1 2 3 4 ones, i.e. the first the number of combinations with 0 ones, etc
  for (i in 1:5) {
    for (j in (i+1):6) {
      for (k in (j+1):7) {
        for (m in (k+1):8) {
          index = (i<5)+(j<5)+(k<5)+(m<5)+1; # represents number of ones in the combination
          numones[index] = numones[index] + 1;
        }
      }
    }
  }
  a = numones/sum(numones);#find number to scale each comb type

  exp4vs4 = log(2.05/1.05)*a[1]*2 + log(1.80/1.30)*a[2]*2;


  expect_equal(t1vs1,exp1vs1, info = "TR001: DSAVEGetTotalVariationFromBulk: 1 vs 1", tolerance=1e-10)
  expect_equal(t4vs4, exp4vs4, info = "TR001: DSAVEGetTotalVariationFromBulk: 4 vs 4", tolerance=1e-10)

})
