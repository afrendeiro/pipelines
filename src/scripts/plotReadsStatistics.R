samples = c(rep("D2_E2F1_R1", 11), rep("D2_E2F7_R1", 11), rep("D2_IN_R1", 11), rep("D2_POL2_R1", 11), rep("D6-F_E2F1_R1", 11), rep("D6-F_E2F7_R1", 11), rep("D6-F_IN_R1", 11), rep("D6-F_POL2_R1", 11), rep("D6-M_E2F1_R1", 11), rep("D6-M_E2F7_R1", 11), rep("D6-M_IN_R1", 11), rep("D6-M_POL2_R1", 11), rep("TB_E2F1_R1", 11), rep("TB_E2F7_R1", 11), rep("TB_IN_R1", 11), rep("TB_POL2_R1", 11))
stages = c(rep("D2", 44), rep("D6-F", 44), rep("D6-M", 44), rep("TB", 44))
types = rep(c("total", "duplicates", "mapped", "paired", "read1", "read2", "proper pair", "itself", "singletons", "diff chr", "diff chr gq"), 16)

path = "/sysdev/s3/share/data/oikopleura/chip-seq/"

S = readLines(paste(path, "flagstat.txt", sep = "/"))

df = data.frame()
for (i in 1:length(S)) {
	df <- rbind(df, data.frame(count = as.numeric(strsplit(x=S[i], split=" ")[[1]][1]), type = types[i], sample = samples[i], stage = stages[i]))
}

# Plot counts
library(ggplot2)
p <- ggplot(data = df, aes(y = count, x = sample)) +
	 geom_bar(stat="identity", width=.5) +
	 facet_grid(type ~ . , scales="free_y") +
	 theme_bw() + theme(strip.text.y = element_text(size = 15))
ggsave(filename=paste(path, "plots/read_counts.png", sep = "/"), plot = p, height = 12, width = 15)

norm <- function(x, i) {
	x / df[df$type == "total","count"][i]
}

# normalize to total (division)
df2 = data.frame()
df2 <- rbind(df2, 
			data.frame(
				count = norm(
							df$count,
							c(rep(1,11),rep(2,11),rep(3,11),rep(4,11),rep(5,11),rep(6,11),rep(7,11),rep(8,11),rep(9,11),rep(10,11),rep(11,11),rep(12,11),rep(13,11),rep(14,11),rep(15,11),rep(16,11))
						),
				type = types,
				sample = samples,
				stage = stages
			)
		)

# round percentages
df2$count <- round(df2$count*100)

# Plot percentages
library(ggplot2)
p <- ggplot(data = df2, aes(y = count, x = sample)) +
	 geom_bar(stat="identity", width=.5) +
	 facet_grid(type ~ .) +
	 ylab("Percentage") +
	 theme_bw() + theme(strip.text.y = element_text(size = 15))
ggsave(filename=paste(path, "plots/read_counts_percentage.png", sep = "/"), plot = p, height = 12, width = 15)