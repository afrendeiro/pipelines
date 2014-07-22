samples = c("D2_E2F1_R1", "D2_E2F7_R1", "D2_IN_R1", "D2_POL2_R1", "D6-F_E2F1_R1", "D6-F_E2F7_R1", "D6-F_IN_R1", "D6-F_POL2_R1", "D6-M_E2F1_R1", "D6-M_E2F7_R1", "D6-M_IN_R1", "D6-M_POL2_R1", "TB_E2F1_R1", "TB_E2F7_R1", "TB_IN_R1", "TB_POL2_R1")
stages = c(rep("D2", 4), rep("D6-F", 4), rep("D6-M", 4), rep("TB", 4))

path = "/sysdev/s3/share/data/oikopleura/chip-seq"

I = read.table(paste(path, "inserts.txt", sep="/"), sep=" ", header = FALSE)
I <- I[seq(1,nrow(I),2)+1,c(4,5)]
colnames(I) <- c("mean", "sd")
I$mean <- as.numeric(gsub(I$mean, pattern = ",", replacement = ""))
I$sd <- as.numeric(gsub(I$sd, pattern = "STD=", replacement = ""))

df = data.frame()
for (i in 1:nrow(I)) {
	df <- rbind(df, data.frame(dist = rnorm(10000, I$mean[i], I$sd[i]), sample = samples[i], stage = stages[i]))
}

library(ggplot2)
p <- ggplot(df, aes(x=dist, fill=sample)) +
	 geom_density(alpha=.3) +
	 facet_grid(stage ~ .) +
	 xlab("Insert size") +
	 theme_bw()

ggsave(filename = paste(path, "plots/insert_size.png", sep="/"), plot = p, height = 5, width = 7)