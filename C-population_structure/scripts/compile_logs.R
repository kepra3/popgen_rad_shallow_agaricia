# Title: Compiling logs and CV error from Admixture results
# Created by: Katharine Prata

# R v4.2.0
library(ggplot2) # v3.4.0

for (vcf_name in c("ac_1diii_nc", "hu_1diii_nc", "lm_1diii_nc-wnr", "all-aga_1diii_nc-wnr")) {
               

  logs <- data.frame(V1 = integer(), V2 = double())
  for (i in 2:10) {
    loglikelihood <- read.table(paste0("../results/admix_runs/20percent_", vcf_name, "_", i, "/loglikelihood.txt"), quote = "\"", comment.char = "")
    subset_log <- loglikelihood[!duplicated(loglikelihood),]
    if (vcf_name != "all-aga_1diii_nc-wnr") {
      subset_log <- subset_log[11:20, ]
    }
    subset_log$V1 <- c(1:10)
    logs <- rbind(logs, subset_log)
  }
  
  logs <- logs[order(logs$V1),]
  
  write.table(logs, file = paste0("../results/admix_runs/logs_", vcf_name, ".txt"), col.names = FALSE, row.names = FALSE, quote = FALSE)
  
  ggplot(logs, aes(x = V1, y = V2)) + geom_line(colour = "red") + geom_point() + theme_bw() +
    ylab("Log-likelihood") +
    scale_x_continuous(name = "K", breaks = c(1:10))  +
    theme(axis.title.x = element_text(size = 12),
          axis.title.y = element_text(size = 12),
          axis.text.x = element_text(size = 8),
          axis.text.y = element_text(size = 8))
  ggsave(paste0("../results/admix_plots/Log-likelihood_", vcf_name, ".pdf"), height = 7, width = 7, units = "cm")
  
  
  CVs <- data.frame(V1 = integer(), V2 = double())
  for (i in 2:10) {
    CV <- read.table(paste0("../results/admix_runs/20percent_", vcf_name, "_", i, "/CV_error.txt"), quote = "\"", comment.char = "")
    subset_CV <- CV[!duplicated(CV),]
    if (vcf_name != "all-aga_1diii_nc-wnr") {
      subset_CV <- subset_CV[11:20, ]
    }
    subset_CV$V1 <- c(1:10)
    CVs <- rbind(CVs, subset_CV)
  }
  
  CVs <- CVs[order(CVs$V1),]
  
  ggplot(CVs, aes(x = V1, y = V4)) + geom_line(colour = "red") + geom_point() + theme_bw() +
    ylab("CV error") +
    scale_x_continuous(name = "K", breaks = c(1:10)) +
    theme(axis.title.x = element_text(size = 12),
          axis.title.y = element_text(size = 12),
          axis.text.x = element_text(size = 8),
          axis.text.y = element_text(size = 8))
  ggsave(paste0("../results/admix_plots/CV_error_", vcf_name, ".pdf"), height = 7, width = 7, units = "cm")
}

