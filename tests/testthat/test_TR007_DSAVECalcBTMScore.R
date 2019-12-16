test_that("TR007 - DSAVECalcBTMScore",{

  extrDir <- downloadData("http://cf.10xgenomics.com/samples/cell-exp/1.1.0/b_cells/b_cells_filtered_gene_bc_matrices.tar.gz", "B10k");
  dataDir = paste0(extrDir,"/filtered_matrices_mex/hg19");
  ds = Read10X(data.dir = dataDir);


  templInfo = DSAVEGetStandardTemplate();

  # Without log transform, sparse matrix. Compare to Matlab results, which are verified using
  # graphs presented in supplementary information of the paper.
  res = DSAVECalcBTMScore(ds, templInfo, silent=TRUE);
  #0.0191 comes from a matlab run; accept some variation
  expect_equal(res$DSAVEScore, 0.0191, info = "TR007: Sparse, not log-transformed", tolerance=0.002)

  #With log transformation:
  #just make sure that the code can run on non-sparse as well
  res2 = DSAVECalcBTMScore(as.matrix(ds[,1:2000]), templInfo, FALSE, 15, TRUE, silent=TRUE);
  #0.0247 comes from a matlab run; accept some variation
  expect_equal(res2$DSAVEScore, 0.0247, info = "TR007: Sparse, not log-transformed", tolerance=0.002)

  })
