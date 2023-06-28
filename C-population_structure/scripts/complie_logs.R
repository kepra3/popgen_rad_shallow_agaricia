library(ggplot2)

# TODO: NOTE only taking half from the log files because accidentally recorded previous results
# If you re-run you need to alter this!

setwd("~/Dropbox/agaricia_project_2019/shalo_ag/gen_project/2 - Population Structure/2b - Admixture/")
for (vcf_name in c("ac_1diii_nc", "hu_1diii_nc", "lm_1diii_nc-wnr")) { #"all-aga_1diii_nc-wnr")) {
               

  logs <- data.frame(V1 = integer(), V2 = double())
  for (i in 2:10) {
    loglikelihood <- read.table(paste0("20percent_", vcf_name, "_", i, "/loglikelihood.txt"), quote = "\"", comment.char = "")
    subset_log <- loglikelihood[!duplicated(loglikelihood),]
    subset_log <- subset_log[11:20, ]
    subset_log$V1 <- c(1:10)
    logs <- rbind(logs, subset_log)
  }
  
  logs <- logs[order(logs$V1),]
  
  write.table(logs, file = paste0("logs_", vcf_name, ".txt"), col.names = FALSE, row.names = FALSE, quote = FALSE)
  
  ggplot(logs, aes(x = V1, y = V2)) + geom_line(colour = "red") + geom_point() + theme_bw() +
    ylab("Log-likelihood") +
    scale_x_continuous(name = "K", breaks = c(1:10))  +
    theme(axis.title.x = element_text(size = 12),
          axis.title.y = element_text(size = 12),
          axis.text.x = element_text(size = 8),
          axis.text.y = element_text(size = 8))
  ggsave(paste0("Log-likelihood_", vcf_name, ".pdf"), height = 7, width = 7, units = "cm")
  
  
  CVs <- data.frame(V1 = integer(), V2 = double())
  for (i in 2:10) {
    CV <- read.table(paste0("20percent_", vcf_name, "_", i, "/CV_error.txt"), quote = "\"", comment.char = "")
    subset_CV <- CV[!duplicated(CV),]
    subset_CV <- subset_CV[11:20, ]
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
  ggsave(paste0("CV_error_", vcf_name, ".pdf"), height = 7, width = 7, units = "cm")
}

