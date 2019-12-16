test_that("TR005 - DSAVEAlignDataset",{
  # First, create a dataset that we know will have a higher UMI distribution
  # for all cells (depends on the current template that comes with the package)
  templInfo = DSAVEGetStandardTemplate();
  extrDir <- downloadData("http://cf.10xgenomics.com/samples/cell-exp/1.1.0/b_cells/b_cells_filtered_gene_bc_matrices.tar.gz", "B10k");
  dataDir = paste0(extrDir,"/filtered_matrices_mex/hg19");
  ds = Read10X(data.dir = dataDir);
  ds2 = ds[, 4001:6000];
  ds3 = DSAVEAlignDataset(ds2, templInfo);
  ds3m = as.matrix(ds3);
  UMIDistr3 = sort(as.vector(colSums(as.matrix(ds3))))
  UMIDistrExp = sort(templInfo$UMIDistr)
  expect_equal(UMIDistr3, UMIDistrExp, info = "TR005: UMIDistr mismatch, no move of counts")

  #Now, create a dataset that we know will not have a higher UMI distr for all cells
  allUMICounts = as.vector(colSums(as.matrix(ds)))
  ds4 = ds[,allUMICounts < 4000];
  ds5 = DSAVEAlignDataset(ds4, templInfo);
  UMIDistr5 = sort(as.vector(colSums(as.matrix(ds5))))
  #first, make sure that the test case has succeeded to create a test where counts
  #need to be moved between cells
  expect_false(isTRUE(all.equal(UMIDistr5, UMIDistrExp)), info = "TR005: UMIDistr matches when it should not, the test case is not testing what it should")
  #then check that it succeeded in the moving, i.e. the sum of all counts should be the same
  expect_equal(sum(UMIDistr5), sum(UMIDistrExp), info = "TR005: UMIDistr false total count")
  #then check that it included the relevant genes only
  expect_equal(sort(rownames(ds5)), sort(templInfo$geneSet), info = "TR005: Failed to synchronize genes")
  #then check that the total number of cells are the same
  expect_equal(dim(ds5)[2], length(templInfo$UMIDistr), info = "TR005: Wrong number of cells")
  #check that sparse input gives sparse output
  expect_true(is(ds4, 'sparseMatrix') & is(ds5, 'sparseMatrix'), info = "TR005: Failure on sparse input")

  #check that non-sparse input gives non-sparse output
  ds6 = as.matrix(ds[,1:2000]);
  ds7 = DSAVEAlignDataset(ds6, templInfo);
  expect_true(is.matrix(ds6) & is.matrix(ds7), info = "TR005: Failure on non-sparse input")


})
