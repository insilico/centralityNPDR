# =====================================================================
# main.R 2024 

# set working directory as needed
setwd("./data")

# https://github.com/insilico/npdr
# library(devtools)
# install_github("insilico/npdr")
library(npdr)
load("df.newz.072421.RData")
load("data.mat.z.072421.RData")
load("class.vec.072421.RData")
load("var.names.072421.RData")

# For easy results later
# Need to read ROIs names
names <- read.csv("roi.names.csv")
# Check:
names[1:5,]

# use the best knn
my.knn <- knnSURF.balanced(df.newz[, 30136], sd.frac = 0.5)

# =====================================================================
# corr diff version of NPDR calls
# =====================================================================
# npdr corr diff lasso
system.time(
  npdr_corr_diff_lasso <- npdr::npdr(outcome=class.vec, dataset=data.mat,
                            regression.type="binomial", attr.diff.type="correlation-data",
                            nbd.method="relieff", nbd.metric = "manhattan", #msurf.sd.frac=0.5,
                            knn=my.knn,
                            use.glmnet = T, glmnet.alpha = 1, #might want to be 0 # set glmnet.alpha=1 for lasso. before 1E-6
                                     #glmnet.lower = -Inf,
                                     glmnet.lam = 0.005, #"lam.1se", #before 0
                                     dopar.nn = F, dopar.reg=F,
                                     corr.attr.names=var.names, padj.method="bonferroni",
                                     verbose = T)
)
dim(npdr_corr_diff_lasso) # 
# Store raw results to file
write.csv(npdr_corr_diff_lasso, "npdr_corr_diff_lasso_2024.csv", row.names=TRUE)

# npdr corr diff non-penalized P-value adjusted
system.time(
  npdr_corr_diff_p_val_adj <- npdr::npdr(outcome=class.vec, dataset=data.mat,
                                         regression.type="binomial", attr.diff.type="correlation-data",
                                         nbd.method="relieff", nbd.metric = "manhattan", #msurf.sd.frac=0.5,
                                         knn=my.knn,
                                         use.glmnet = F, glmnet.alpha = 1E-6, #might want to be 0
                                         glmnet.lower = -Inf,
                                         glmnet.lam = 0,
                                         dopar.nn = F, dopar.reg=F,
                                         corr.attr.names=var.names, padj.method="bonferroni",
                                         verbose = T)
)
dim(npdr_corr_diff_p_val_adj) # 
# Store raw results to file
write.csv(npdr_corr_diff_p_val_adj, "npdr_corr_diff_p_val_adj_2024.csv", row.names=FALSE)

# remove sex and (mdd) status from data
df.newz.rmsex <- df.newz[,-30137]
df.newz.rmsex.rmstatus <- df.newz.rmsex[,-30136]
# =====================================================================
# numeric abs version of NPDR calls
# =====================================================================
# npdr numeric abs lasso
system.time(
  npdr_numeric_abs_lasso <- npdr::npdr(outcome=class.vec, dataset=df.newz.rmsex.rmstatus,
                           regression.type="binomial", attr.diff.type="numeric-abs",
                           nbd.method="relieff", nbd.metric = "manhattan", #msurf.sd.frac=0.5,
                                     knn=my.knn,
                                     use.glmnet = T, glmnet.alpha = 1, #might want to be 0 # set glmnet.alpha=1 for lasso. before 1E-6
                                     #glmnet.lower = -Inf,
                                     glmnet.lam = 0.005, #"lam.1se", #before 0
                                     dopar.nn = F, dopar.reg=F,
                                     corr.attr.names=var.names, padj.method="bonferroni",
                                     verbose = T)
)

npdr_numeric_abs_lasso[1:10,1] # 
dplyr::slice_max(npdr_numeric_abs_lasso,order_by=scores,n=20)

# Store raw results to file
write.csv(npdr_numeric_abs_lasso, "npdr_numeric_abs_lasso_2024.csv", row.names=TRUE)

# npdr numeric abs ridge
system.time(
  npdr_numeric_abs_ridge <- npdr::npdr(outcome=class.vec, dataset=df.newz.rmsex.rmstatus,
                            regression.type="binomial", attr.diff.type="numeric-abs",
                            nbd.method="relieff", nbd.metric = "manhattan", #msurf.sd.frac=0.5,
                                       knn=my.knn,
                                       use.glmnet = T, glmnet.alpha = 1E-6, #might want to be 0
                                       glmnet.lower = -Inf,
                                       glmnet.lam = .005, # maybe try something not 0 #lam_1se
                                       dopar.nn = F, dopar.reg=F,
                                       corr.attr.names=var.names, padj.method="bonferroni",
                                       verbose = T)
)
dim(npdr_numeric_abs_ridge) # 
# Store raw results to file
write.csv(npdr_numeric_abs_ridge, "npdr_numeric_abs_ridge_2024.csv", row.names=FALSE)


# npdr numeric abs non-penalized P-value adjusted
system.time(
  npdr_numeric_abs_p_val_adj <- npdr(df.newz[,30136],df.newz[,1:30135], 
                                     regression.type="binomial", attr.diff.type="numeric-abs",
                                     nbd.method="relieff", nbd.metric = "manhattan", 
                                     knn=my.knn, msurf.sd.frac=.5, 
                                     use.glmnet = F, dopar.nn = F, dopar.reg = F, 
                                     padj.method="bonferroni", verbose=T)
)
dim(npdr_numeric_abs_p_val_adj) 
# Store raw results to file
write.csv(npdr_numeric_abs_p_val_adj, "npdr_numeric_abs_p_val_adj_2024.csv", row.names=FALSE)

# =====================================================================
# utiliy function
# A complete function that converts the results of "numeric-abs"
# version of npdr into complete names results, with pvalue
# This is necessary to create graphs and then apply degree centrality
# =====================================================================
t5h.full.results <- function (res, name_col, pval_col) {
  # name_col is the number of column that stores attr name
  # pval_col is the number of column that stores pvalue
  library(stringr)
  strs <- res[, name_col]
  a <- substring(strs, 18, last=1000000L)
  b <- str_split(a, '_')
  first <- sapply(b, '[[', 1)
  second <- sapply(b, '[[', 2)
  alist <- list(first, second)
  mat <- do.call(cbind, alist)
  full.res <- cbind(as.numeric(mat[, 1]), as.numeric(mat[, 2]))
  dim(full.res)
  cn <- names[as.numeric(full.res[, 1]), 2]
  gn <- names[as.numeric(full.res[, 1]), 3]
  full.res <- cbind(full.res, cn, gn)
  
  cn <- names[as.numeric(full.res[, 2]), 2]
  gn <- names[as.numeric(full.res[, 2]), 3]
  full.res <- cbind(full.res, cn, gn)
  
  full.res <- cbind(full.res, res[, pval_col])	
  
  return(full.res)
}

# =====================================================================
# npdr_numeric_abs_lasso then applying degree centrality
# =====================================================================
npdr_numeric_abs_lasso[1:10,]
npdr_numeric_abs_lasso$scores

# change column name
names(npdr_numeric_abs_lasso)[1] <- "roipairs"
lasso_vars <- npdr_numeric_abs_lasso
lasso_vars_sorted <- sort(lasso_vars, decreasing = T)
lasso_vars_sorted[1:10,]
# =====================================================================
# Now try degree centrality
#### rf centrality methods####
roi.pairs.str <- lasso_vars_sorted[1:100]
roi.pairs.edges <- do.call(rbind,names(roi.pairs.str)) 

# turn into graph 
roi_graph <- graph_from_edgelist(roi.pairs.edges, directed=F) 
plot(roi_graph)
#### degree ####
numeric_abs_lasso_degree <- sort(degree(roi_graph),decreasing=T)  
numeric_abs_lasso_degree
write.csv(numeric_abs_lasso_degree, "numeric_abs_lasso_deg_2024.csv", 
          row.names=FALSE)



# library(igraph)
grph <- graph_from_edgelist(nums, directed=F)
plot(grph)
grph
length(grph)
str(grph[1])
names(grph[1])
names(grph[2])
grph[1]
grph[2]

# Compute the degree centrality
deg.cen <- centr_degree(grph, "all")
str(deg.cen)
deg.cen$res

# An attempt to match the name to the position
len <- length(deg.cen$res)
degcen <- matrix("", len, 2)
for (i in 1:len) {
  degcen[i, 1] <- names(grph[1])[i]
  degcen[i, 2] <- deg.cen$res[i]
}
degcen[1:10, ]

# Sort on the second column
degcen.sorted <- degcen[order(as.numeric(degcen[, 2]), decreasing=TRUE),]
degcen.sorted[1:20, ]

degcen.sorted[1:20, 1]

# Store degree centrality for np cov 2k
write.csv(degcen.sorted, "npdr_numeric_abs_lasso_deg_2024.csv", 
          row.names=FALSE)

# =====================================================================
# npdr_numeric_abs_p_val_adj then applying degree centrality
# =====================================================================
# filter to take the top 200 beta.z.attr because pairs, turn that into the network
converted_npdr_numeric_abs_p_val_adj <- t5h.full.results(npdr_numeric_abs_p_val_adj, 1, 2)

# ask McK for help
# Prepare the ROI graph
# Obtain input in the form of matrix 1011 by 2 (pairs of indices)
#np.2num <- converted_npdr_numeric_abs_p_val_adj[1:1011, 1:2]
#np.2num <- converted_npdr_numeric_abs_p_val_adj[1:1011, c(3,5)]
np.2num <- converted_npdr_numeric_abs_p_val_adj[1:200, c(3,5)]

# Change name to make it the same as in the sample code for the method
nums <- np.2num

# Checking display
#nums[1:5,]

library(igraph)
roi_graph.t5h <- graph_from_edgelist(nums, directed=F)
plot(roi_graph.t5h)
library(tcltk)

# color the nodes by louvain clustering
lou <- cluster_louvain(roi_graph.t5h)
table(lou$membership)  # has more clusters than greedy
tkplot(roi_graph.t5h, 
       vertex.color=lou$membership ,      
       vertex.label.dist=0.8, vertex.label.color="blue",
       vertex.label.cex=0.5,
       vertex.size=degree(roi_graph.t5h)  # size by main effect
)

# better graph
tkplot(roi_graph.t5h, 
       vertex.color="red",      
       vertex.label.dist=1.5, vertex.label.color="blue",
       vertex.size=degree(roi_graph.t5h)  # size by main effect
)
#save converted
write.csv(converted_npdr_numeric_abs_p_val_adj, "converted_npdr_numeric_abs_p_val_adj_april_2024.csv", 
          row.names=FALSE)
#save original jic 
write.csv(npdr_numeric_abs_p_val_adj, "npdr_numeric_abs_p_val_adj_april_2024.csv", 
          row.names=FALSE)
# Try coloring the nodes by louvain clustering
lou <- cluster_louvain(roi_graph.t5h)
table(lou$membership)  # has more clusters than greedy
tkplot(roi_graph.t5h, 
       vertex.color=lou$membership ,      
       vertex.label.dist=1.5, vertex.label.color="blue",
       vertex.size=degree(roi_graph.t5h)  # size by main effect
)
# change font size
lou <- cluster_louvain(roi_graph.t5h)
table(lou$membership)  # has more clusters than greedy
tkplot(roi_graph.t5h, # .cex = 0.09 then 0.9
       vertex.color=lou$membership ,      
       vertex.label.dist=1.0, vertex.label.color="blue",vertex.label.cex=0.75,
       vertex.size=degree(roi_graph.t5h)  # size by main effect
)
# change font size bigger
lou <- cluster_louvain(roi_graph.t5h)
table(lou$membership)  # has more clusters than greedy
tkplot(roi_graph.t5h, # .cex = 0.09 then 0.9
       vertex.color=lou$membership ,      
       vertex.label.dist=1.0, vertex.label.color="blue",vertex.label.cex=0.85,
       vertex.size=degree(roi_graph.t5h)  # size by main effect
)
# change font size bigger
lou <- cluster_louvain(roi_graph.t5h)
table(lou$membership)  # has more clusters than greedy
tkplot(roi_graph.t5h, # .cex = 0.09 then 0.9
       vertex.color=lou$membership ,      
       vertex.label.dist=1.0, vertex.label.color="blue",vertex.label.cex=1.5,
       vertex.size=degree(roi_graph.t5h)  # size by main effect
)

# =====================================================================
# Now try degree centrality
# library(igraph)
grph <- graph_from_edgelist(nums, directed=F)
plot(grph)
grph
length(grph)
str(grph[1])
names(grph[1])
names(grph[2])
grph[1]
grph[2]

# Compute the degree centrality
deg.cen <- centr_degree(grph, "all")
str(deg.cen)
deg.cen$res

# An attempt to match the name to the position
len <- length(deg.cen$res)
degcen <- matrix("", len, 2)
for (i in 1:len) {
  degcen[i, 1] <- names(grph[1])[i]
  degcen[i, 2] <- deg.cen$res[i]
}
degcen[1:10, ]

# Sort on the second column
degcen.sorted <- degcen[order(as.numeric(degcen[, 2]), decreasing=TRUE),]
degcen.sorted[1:20, ]

degcen.sorted[1:20, 1]

# Store degree centrality for np cov 2k
write.csv(degcen.sorted, "npdr_numeric_abs_p_val_adj_deg_2024.csv", 
          row.names=FALSE)

# npdr_numeric_abs_p_val_adj # take top 200

# Change the results to the 2names format:
# (reg = diagnosis as outcome)
converted_npdr_numeric_abs_p_val_adj <- t5h.full.results(npdr_numeric_abs_p_val_adj, 1, 2)

# Prepare the ROI graph
# Obtain input in the form of matrix 1011 by 2 (pairs of indices)
np.2num <- converted_npdr_numeric_abs_p_val_adj[1:1011, 1:2]

# Change name to make it the same as in the sample code for the method
nums <- np.2num
#dim(nums)

# Checking display
#nums[1:5,]

roi_graph.t5h <- graph_from_edgelist(nums, directed=F)
plot(roi_graph.t5h) 

# degree centrality setup
grph <- graph_from_edgelist(nums, directed=F)
plot(grph)
grph
length(grph)
str(grph[1])
names(grph[1])
names(grph[2])
grph[1]
grph[2]

# Compute the degree centrality
deg.cen <- centr_degree(grph, "all")
str(deg.cen)
deg.cen$res

# An attempt to match the name to the position
len <- length(deg.cen$res)
degcen <- matrix("", len, 2)
for (i in 1:len) {
  degcen[i, 1] <- names(grph[1])[i]
  degcen[i, 2] <- deg.cen$res[i]
}
degcen[1:10, ]

# Sort on the second column
degcen.sorted <- degcen[order(as.numeric(degcen[, 2]), decreasing=TRUE),]
degcen.sorted[1:20, ]

degcen.sorted[1:20, 1]

# Store degree centrality for np cov 2k
write.csv(degcen.sorted, "deg.npdr_numeric_abs_p_val_adj.2024.csv", 
          row.names=FALSE)

# =====================================================================
# npdr_numeric_abs_ridge then applying degree centrality
# =====================================================================
# npdr_numeric_abs_ridge # take top 200

# Prepare the ROI graph
# Obtain input in the form of matrix 1011 by 2 (pairs of indices)
np.2num <- npdr_numeric_abs_ridge[1:1011, 1:2]

# Change name to make it the same as in the sample code for the method
nums <- np.2num
#dim(nums)

# Checking display
#nums[1:5,]

library(igraph)
roi_graph.t5h <- graph_from_edgelist(nums, directed=F)
plot(roi_graph.t5h) 

# degree centrality setup
grph <- graph_from_edgelist(nums, directed=F)
plot(grph)
grph
length(grph)
str(grph[1])
names(grph[1])
names(grph[2])
grph[1]
grph[2]

# Compute the degree centrality
deg.cen <- centr_degree(grph, "all")
str(deg.cen)
deg.cen$res

# An attempt to match the name to the position
len <- length(deg.cen$res)
degcen <- matrix("", len, 2)
for (i in 1:len) {
  degcen[i, 1] <- names(grph[1])[i]
  degcen[i, 2] <- deg.cen$res[i]
}
degcen[1:10, ]

# Sort on the second column
degcen.sorted <- degcen[order(as.numeric(degcen[, 2]), decreasing=TRUE),]
degcen.sorted[1:20, ]

degcen.sorted[1:20, 1]

# Store degree centrality for np cov 2k
write.csv(degcen.sorted, "deg.npdr_numeric_abs_ridge.2024.csv", 
          row.names=FALSE)

# =====================================================================
# Range R Random forest
# =====================================================================
# install.packages("ranger")
library(ranger)

df.newz.rmsex <- df.newz[,-30137]
#rm(df.newz)
df.newz.rmsex$status <- as.factor(df.newz.rmsex$status)
rf_results <- ranger(x=df.newz.rmsex[,-30136], y=df.newz.rmsex[,30136],
                     importance = "permutation", num.trees=5000,
                     classification = T)

rf_vars <- rf_results$variable.importance
rf_vars_sorted <- sort(rf_vars, decreasing = T)
rf_vars_sorted[1:10]
1-rf_results$prediction.error

#write.csv(rf_results$variable.importance, "random_forest_2024.csv", 
#row.names=FALSE)

# ask McK for help
#### rf centrality methods####
roi.pairs.str <- rf_vars_sorted[1:200]
roi.pairs.edges <- do.call(rbind,strsplit(names(roi.pairs.str)[1:200],"[_]")) 
# remove unecessary columns
roi.pairs.edges <- roi.pairs.edges[,-1]
roi.pairs.edges <- roi.pairs.edges[,-1]
roi.pairs.edges <- roi.pairs.edges[,-1]
roi.pairs.edges
# turn into graph 
roi_graph <- graph_from_edgelist(roi.pairs.edges, directed=F) 
plot(roi_graph)

# size nodes by degree centrality
plot(roi_graph,vertex.size=(degree(roi_graph)),
     vertex.color="red", label.color="blue", vertex.label.cex=0.9)

# better graph
tkplot(roi_graph, 
       vertex.color="red",      
       vertex.label.dist=1.5, vertex.label.color="blue",
       vertex.size=degree(roi_graph)  # size by main effect
)
# Try coloring the nodes by louvain clustering
lou <- cluster_louvain(roi_graph)
table(lou$membership)  # has more clusters than greedy
tkplot(roi_graph, 
       vertex.color=lou$membership ,      
       vertex.label.dist=1.5, vertex.label.color="blue",
       vertex.size=degree(roi_graph)  # size by main effect
)

#### degree ####
rf_degree <- sort(degree(roi_graph),decreasing=T)  
rf_degree
write.csv(rf_degree, "random_forest_deg_2024.csv", 
          row.names=FALSE)

# get indices  
rf_degree_idx <- do.call(rbind,strsplit(names(rf_degree),"ROI"))[,2] 
names(rf_degree)
#rf_degree_idx <- do.call(rbind,names(rf_degree),"ROI")[,2]

#scores
scores_all = as.numeric(rf_degree) 

#variable names of the scores
varnames_all = rf_degree_idx 

#not sure if this works
x = rbind(scores_all,varnames_all)

write.csv(x, "random_forest_deg_2024.csv", 
          row.names=FALSE)

################## variation inflation factor (vif) #############
#install.packages("spatialRF")
library(spatialRF)
rf_vars <- rf_results$variable.importance
rf_vars_sorted <- sort(rf_vars, decreasing = T)
rf_vars_sorted[1:10]
which_rfimp_top <- which(rf_results$variable.importance>=rf_vars_sorted[200])
temp_rf <- df.newz.rmsex.rmstatus[,which_rfimp_top]
rf_vif <- vif(temp_rf)
rf_vif_rois <- t5h.full.results(rf_vif, 1, 2) 
rf_vif_rois[,c(3,5,7)]

# npdr importances are pre-sorted
which_npdr_top <- which(colnames(df.newz.rmsex.rmstatus) %in% rownames(npdr_numeric_abs_lasso)[2:201]) # skip intercept
temp_npdr <- df.newz.rmsex.rmstatus[,which_npdr_top]
npdr_vif <- vif(temp_npdr)
npdr_vif_rois <- t5h.full.results(npdr_vif, 1, 2) 
npdr_vif_rois[,c(3,5,7)]



##### AAl Atlas (90 ROIs) #####

# delete
library(ranger)
# Henry's call
rf_results <- ranger(class.vec~., data = data.upper.mat, # need to adjust these args
                     importance = "permutation", num.trees=5000,
                     classification = T)
# my previous randomForest call:
# rf_results <- randomForest(df.newz[,1:30135], as.factor(df.newz[,30136]), ntree=5000) 
## my guess, gives an error
rf_results <- ranger((as.factor(df.newz[,30136]))~., data=df.newz[,1:30135], 
                     importance = "permutation", num.trees=100,
                     classification = T)
# Error in formula.default(formula) : invalid formula
## my second guess, gives an error
rf_results <- ranger((as.factor(df.newz[,30136]))~., data=data.mat, 
                     importance = "permutation", num.trees=5000,
                     classification = T)

rf_results <- ranger(class.vec~., data=data.mat, 
                     importance = "permutation", num.trees=500,
                     classification = T)

# Error in terms.formula(f, data = data) : 
# duplicated name 'SFG_L_7_1.SFG_R_7_1' in data frame using '.'

rf_vars <- rf_results$variable.importance 

rf_vars_sorted <- sort(rf_vars, decreasing = T)

#### rf centrality methods####
roi.pairs.str <- rf_vars_sorted[1:100]
roi.pairs.edges <- do.call(rbind,strsplit(names(roi.pairs.str)[1:100],"[.]")) 

# turn into graph 
roi_graph <- graph_from_edgelist(roi.pairs.edges, directed=F) 

#### degree ####
npdr_degree <- sort(degree(roi_graph),decreasing=T)  

# =====================================================================
# Random Forest (Liz's using RandomForest)
# =====================================================================
library(randomForest)
library(dplyr)
library(tibble)

set.seed(1210)

# Run the random forest method for the correlations between ROIs. 
# We use the non-default values of ntree

# First execution
t500.beg <- Sys.time()
t500.fit.1210 <- randomForest(df.newz[,1:30135], as.factor(df.newz[,30136]), 
ntree=5000) 
t500.end <- Sys.time()
t500.time <- t500.end - t500.beg
t500.time

# Print the object returned by running the method
# print(t500.fit.1210)

# Get the most important ROIs
t500.importance.1210 <- data.frame(importance(t500.fit.1210))
t500.results.1210 <- t500.importance.1210 %>% rownames_to_column() %>% 
  arrange(-MeanDecreaseGini)

# Get the current 200 of most important pairs
ranfor.t500.top200.1210 <- t500.results.1210[1:200,]
dim(ranfor.t500.top200.1210)

ranfor.t500.top200.1210[1:20,1]

write.table(ranfor.t500.top200.1210[,1], "random_forest_top200_2024.txt", 
            row.names=FALSE, col.names=FALSE, sep="\n")

library(randomForest)
library(dplyr)
library(tibble)

set.seed(1210)

# Run the random forest method for the correlations between ROIs. 
# Using non-default value ntree = 5000 

# First execution
btime <- Sys.time()
rf.1210 <- randomForest(df.newz[,1:30135], as.factor(df.newz[,30136]), 
                        ntree=5000) 
etime <- Sys.time()
tim <- etime - btime
tim

# Get the most important ROIs
rf.imp.1210 <- data.frame(importance(rf.1210))
rf.results.1210 <- rf.imp.1210 %>% rownames_to_column() %>% 
  arrange(-MeanDecreaseGini) # use permutation

# Get the current 200 of most important pairs 
rf.top2k.1210 <- rf.results.1210[1:200,]
dim(rf.top2k.1210)

rf.top2k.1210[1:10,1]

write.csv(rf.top2k.1210, "rf_top200_2024.csv", row.names=F)

# Converting into 'pairs' format
rf.t2k.num <- substring(rf.top2k.1210[,1], 18, last=1000000L)
length(rf.t2k.num)
rf.t2k.num[1:10]

library(stringr)

# Fragment of code that converts the "xxx_yyy" format into "xxx" "yyy"
# Now get just the ROIs numbers
rf.t2k.roi.nums <- str_split(rf.t2k.num, "_")
rf.t2k.roi.nums[1:20]
length(rf.t2k.roi.nums)

rf.t2k.first.num <- sapply(rf.t2k.roi.nums, "[[", 1)
rf.t2k.second.num <- sapply(rf.t2k.roi.nums, "[[", 2)

rf.t2k.first.num[1:20]
rf.t2k.second.num[1:20]
# This is good

# Make an intermediary vector
rf.t2k.roi.nums[1:20]
rf.t2k.roi.vnums <- unlist(rf.t2k.roi.nums)
rf.t2k.roi.vnums[1:20]
length(rf.t2k.roi.vnums)

# Finally make it into a matrix of numbers
rf.t2k.roi.pairs <- matrix(0, nrow=200, ncol=2)
for (i in 1:200) {
  rf.t2k.roi.pairs[i, 1] = rf.t2k.roi.vnums[2*i - 1]
  rf.t2k.roi.pairs[i, 2] = rf.t2k.roi.vnums[2*i]
}
dim(rf.t2k.roi.pairs)
# Store to file
write.csv(rf.t2k.roi.pairs, "rf.1210.t2k.5k.pairs.2024.csv", 
          row.names=FALSE)

# Convert to a single number
rf.5k <- matrix(0, 200, 1)
for (i in 1:200) {
  rf.5k[i] <- 246 * (as.numeric(rf.t2k.roi.pairs[i, 1]) - 1) +
    as.numeric(rf.t2k.roi.pairs[i, 2])
}

rf.5k[1:20]

write.csv(rf.5k, "rf.5k.1210.sglenum.2024.csv", row.names=F)
