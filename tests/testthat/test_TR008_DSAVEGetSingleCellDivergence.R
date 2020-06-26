test_that("TR008 - DSAVEGetSingleCellDivergence",{

  d <- rbind(c(0,1,2,1),c(3,3,1,2));
  ds <- as.matrix(d);

  divs = DSAVEGetSingleCellDivergence(ds, 3, silent = TRUE)$divs;

  ds11 <- matrix(1,30000,2);
  divs11 = DSAVEGetSingleCellDivergence(ds11, 400, silent = TRUE)$divs;


  d2 <- rbind(c(0,14,2,1),c(17,0,1,2));
  ds2 <- as.matrix(d2);

  divs2 = DSAVEGetSingleCellDivergence(ds2, 3, silent = TRUE)$divs;

  d3 <- rbind(c(0,3,2,1),c(3,0,1,2));
  ds3 <- as.matrix(d3);

  divs3 = DSAVEGetSingleCellDivergence(ds3, 3, silent = TRUE)$divs;

  d4 <- rbind(c(0,14,5,1),c(17,0,1,2));
  ds4 <- as.matrix(d4);
  divs4 = DSAVEGetSingleCellDivergence(ds4, 6, silent = TRUE)$divs;

  #check that more divergent cells get higher divergence
  expect_true((divs[4] < divs[3]) & (divs[4] < divs[1]) &&
              (divs[4] <= divs[2]), info = "TR008: Value not ok")

  #Also explicitly check the values against matlab.
  expVals = c(1.1240803,0.8139254,1.6023828,0.8139254);
  expect_equal(divs, expVals, info = "TR008: Values not same as in matlab", tolerance=1e-4)

  #check that downsampling gives equal results for the case set up
  expect_equal(divs2, divs3, info = "TR008: Down-sampling not ok", tolerance=1e-4)

  #check that NA is given for only the cells with too few counts
  expect_true(is.na(divs4[4]) & !is.na(divs4[3]), info = "TR008: NA not ok")

  #also make sure it runs for a real dataset
  extrDir <- downloadData("http://cf.10xgenomics.com/samples/cell-exp/1.1.0/b_cells/b_cells_filtered_gene_bc_matrices.tar.gz", "B10k");
  dataDir = paste0(extrDir,"/filtered_matrices_mex/hg19");
  ds10 = Read10X(data.dir = dataDir);

  ds10 = ds10[, 1:100];
  divs10 = DSAVEGetSingleCellDivergence(ds10, 200, silent = TRUE)$divs;
  #check that we have no NA
  expect_equal(sum(is.na(divs10)), 0, info = "TR008: NA Values in b10k")

  #Check the gene-wise divergence
  dgs <- cbind(c(8,0),c(0,8), c(8,0), c(4,0));
  dsgs <- as.matrix(dgs);
  row.names(dsgs) = c("Gene A", "Gene B")
  colnames(dsgs) = c("Cell A", "Cell B", "Cell C", "Cell D")

  logFacs = GetLogFactorials(2000);
  rs = rowMeans(dsgs)
  prob = c(0.75,0.25)
  a = LogBinomialPDF(dsgs[,2], prob = prob, logFacs);
  expRes1 = dbinom(8, 8, 0.25, log = TRUE)
  expect_equal(a[1], expRes1, info = "TR008: LogBinomialPDF not ok", tolerance=1e-4, check.names=F)

  geneDivs = DSAVEGetSingleCellDivergence(dsgs, 4, silent = TRUE)$geneDivs;
  expRes = -dbinom(4, 4, 0.25, log = TRUE)
  expect_equal(geneDivs[2,2], expRes, info = "TR008: Gene-wise divergence not ok", tolerance=1e-4, check.names=F)
  })





