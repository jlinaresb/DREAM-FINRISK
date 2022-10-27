require(data.table)
require(ggplot2)
require(ggpubr)
require(viridis)

setwd(file.path(here::here(), "src/jlb_m1/results/"))
files <- list.files()

experiment_id <- "pca"
experiment <- readRDS(files[grep(experiment_id, files)])

pred <- experiment$prediction
times <- seq(1, 12, 1)
methods <- names(experiment$models)
subset <- "intval" # train intval extval

s <- ifelse(subset == "train", 1,
        ifelse(subset == "intval", 2,
          ifelse(subset == "extval", 3, NA)))

c <- list()
b <- list()
d <- list()

for (i in seq_along(pred)) {
  # Get C-Index
  c[[i]] <- as.data.frame(pred[[i]][s])[1, ]
  # Get Brier error
  b[[i]] <- as.data.frame(pred[[i]][s])[2, ]
  # Get Delta-Index
  d[[i]] <- as.data.frame(pred[[i]][s])[3, ]
}

res <- list()
res$CIndex <- as.data.frame(rbindlist(c))
rownames(res$CIndex) <- times
colnames(res$CIndex) <- methods

res$Briererror <- as.data.frame(rbindlist(b))
rownames(res$Briererror) <- times
colnames(res$Briererror) <- methods

res$DeltaCIndex <- as.data.frame(rbindlist(d))
rownames(res$DeltaCIndex) <- times
colnames(res$DeltaCIndex) <- methods

predictions <- list(CIndex = res$CIndex,
                   Briererror = res$Briererror,
                   DeltaCIndex = res$DeltaCIndex)

# C-Index
data <- data.table::melt(predictions$CIndex)
data$times <- rep(times, length(methods))
p1 <- ggplot(data = data, aes(x = times, y = value, group = variable)) +
  geom_line(aes(color = variable)) + geom_point(aes(color = variable)) +
  scale_color_manual(values = viridis(15)) + ylim(0, 1) +
  theme_bw() +
  ggtitle(paste0("C-Index - ", subset))

p1
# Brier error
data <- data.table::melt(predictions$Briererror)
data$times <- rep(times, length(methods))
p2 <- ggplot(data = data, aes(x = times, y = value, group = variable)) +
  geom_line(aes(color = variable)) + geom_point(aes(color = variable)) +
  scale_color_manual(values = viridis(15)) + ylim(0, 1) +
  theme_bw() +
  ggtitle(paste0("Brier error - ", subset))
p2

# Delta Index
data <- data.table::melt(predictions$DeltaCIndex)
data$times <- rep(times, length(methods))
p3 <- ggplot(data = data, aes(x = times, y = value, group = variable)) +
  geom_line(aes(color = variable)) + geom_point(aes(color = variable)) +
  scale_color_manual(values = viridis(15)) + ylim(0, 1) +
  theme_bw() +
  ggtitle(paste0("Delta Index - ", subset))
p3

plot <- ggarrange(p1, p2, p3, nrow = 3, ncol = 1, common.legend = TRUE)
