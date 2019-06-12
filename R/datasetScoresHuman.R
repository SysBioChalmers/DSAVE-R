#' datasetScoresHuman
#'
#' DSAVE scores from public datasets evaluated with the standard template.
#' Can be used for comparison against the score generated to see if it is
#' comparable to scores from other datasets.
#'
#' @format list of lists, where the inner list contains the dataset name (as named
#' in the DSAVE publication) and the DSAVE score.
#' @source Calculated with the DSAVE-M package in MATLAB, but the results should
#' be the same.
#' \tabular{lll}{
#' 	BC \tab Azizi, E. et al. Single-Cell Map of Diverse Immune Phenotypes in the Breast Tumor Microenvironment. Cell (2018). doi:10.1016/j.cell.2018.05.060. \tab 50/50 mix of T cells from the two samples 'BC1_BLOOD' and 'BC4_BLOOD' from the BC dataset, in total 3,532 cells. \cr
#' 	OC \tab	Schelker, M. et al. Estimation of immune cell content in tumour tissue using single-cell RNA-seq data. Nat. Commun. 8, 2032 (2017). \tab All macrophages from the OC dataset, in total 2,057 cells. \cr
#'	LC \tab	Lambrechts, D. et al. Phenotype molding of stromal cells in the lung tumor microenvironment. Nat. Med. 24, 1277–1289 (2018). \tab All t cells from the healthy tissue of patient 3, 4 and 5 from the LC dataset, in total 3,471 cells. \cr
#'	LIVC \tab Zheng, C. et al. Landscape of Infiltrating T Cells in Liver Cancer Revealed by Single-Cell Sequencing. Cell 169, 1342-1356.e16 (2017). \tab All cells from the LIVC dataset, in total 4,070 cells. Note that UMIs are not used for this dataset, which could explain the higher variation. \cr
#'	PBMC68k \tab	Zheng, G. X. Y. et al. Massively parallel digital transcriptional profiling of single cells. Nat. Commun. 8, 14049 (2017). \tab ~68,000 T cells. All T cells from the PBMC68k dataset, in total 48,657 cells. \cr
#'	B10k \tab Zheng, G. X. Y. et al. Massively parallel digital transcriptional profiling of single cells. Nat. Commun. 8, 14049 (2017). \tab ~10,000 FACS-sorted B cells (CD19+). All cells from the B10k dataset, in total 10,085 cells. \cr
#'	CD4TMEM \tab	Zheng, G. X. Y. et al. Massively parallel digital transcriptional profiling of single cells. Nat. Commun. 8, 14049 (2017). \tab ~10,000 FACS-sorted CD4+ T memory cells. All cells from the CD4TMEM dataset, in total 10,224 cells. \cr
#'	HCA CB \tab Li, B. et al. Census of Immune Cells. Human Cell Atlas Data Portal (2018). \tab All T cells from patient 'CB1', 'CB2' and 'CB3' from the HCA CB dataset, in total 30,462 cells. \cr
#'	CD8T \tab Chen, J. et al. PBMC fixation and processing for Chromium single-cell RNA sequencing. J. Transl. Med. 16, 198 (2018). \tab All cells from the TCD8 dataset, in total 5,662 cells. \cr
#' }
"datasetScoresHuman"
