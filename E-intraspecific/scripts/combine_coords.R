
setwd("~/Dropbox/agaricia_project_2019/shalo_ag/gen_project/3 - Spatial/3b - Within & between plots/spatial-agaricia/")

WP05 <- read.csv("data/scaled_HORIZ_annotations_cur_kal_05m_20200214.csv")
WP10 <- read.csv("data/scaled_HORIZ_annotations_cur_kal_10m_20200214.csv")
WP20 <- read.csv("data/scaled_HORIZ_annotations_cur_kal_20m_20200214.csv")
CA05 <- read.csv("data/scaled_HORIZ_annotations_cur_cas_05m_20201212.csv")
CA10 <- read.csv("data/scaled_HORIZ_annotations_cur_cas_10m_20201212.csv")
CA20 <- read.csv("data/scaled_HORIZ_annotations_cur_cas_20m_20201212.csv")
SB05 <- read.csv("data/scaled_HORIZ_annotations_cur_sna_05m_20200303.csv")
SB10 <- read.csv("data/scaled_HORIZ_annotations_cur_sna_10m_20200303.csv")
SB20 <- read.csv("data/scaled_HORIZ_annotations_cur_sna_20m_20200303.csv")
SQ12 <- read.csv("data/scaled_HORIZ_annotations_cur_seb_10m_20201210.csv")
SQ20 <- read.csv("data/scaled_HORIZ_annotations_cur_seb_20m_20201210.csv")

all_plots <- rbind(WP05, WP10, WP20, CA05, CA10, CA20, SB05, SB10, SB20, SQ12, SQ20)
colnames(all_plots) <- c("Individual", "x", "y", "z", "range")
write.csv(all_plots, file = "data/coord_horiz_metdata.csv", quote = FALSE, row.names = FALSE)
