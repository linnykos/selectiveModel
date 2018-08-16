dat <- read.csv("../experiments/profile_results.txt", sep = ",")

rownames(dat) <- dat[,1]
dat <- dat[,-1]
colnames(dat) <- as.character(c(4,6,9,12,16,24))

total_remainder <- dat["total",] - dat["interval_list",]
interval_remainder <- dat["interval_list",] - colSums(dat[grep("\\.", rownames(dat)),])

dat <- rbind(total_remainder, interval_remainder, dat[-c(1,2),])

col_vec <- c("gray", "gold", "coral4", "coral1", "coral3", "coral", "chocolate4",
             "coral2", "blue")

dat <- as.table(as.matrix(dat))
barplot(dat, col = col_vec)
