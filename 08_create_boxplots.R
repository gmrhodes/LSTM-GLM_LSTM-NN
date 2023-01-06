# This file creates the boxplots of losses and c-indexes for the simulation.

############################### USER SET-UP ##############################
#Survival model ("aft" or "cox")
surv_model = "cox"

#Directory path of folder containing results from 07_loss_cIndex.R
##Must contain sub-directories "aft" and "cox"
result_path = "D:/Simulations/results/"


############################### CREATE LOSS BOXPLOTS ##############################
#Read in training & testing losses
trainLoss_mat = read.table(paste0(result_path, surv_model, "/train_loss.csv"), sep=",", header=T)
testLoss_mat = read.table(paste0(result_path, surv_model, "/test_loss.csv"), sep=",", header=T)

#Y-axis upper-limit
yUpper = ifelse(surv_model=="aft", 6, 6.5)

#Plot training losses
png(paste0(result_path, surv_model, "/trainLoss_plot.png"))
med = apply(trainLoss_mat, 2, median)
colfunc = colorRampPalette(c("chartreuse4","chartreuse2","yellow","orange","red"))
colVec = colfunc(8)[rank(med)]
par(cex.axis=1.25, cex.main=1.5) 
boxplot(as.list(as.data.frame(trainLoss_mat)), main=paste0(toupper(surv_model), ": Training Loss"), ylab="",
        names=c("B", "L", "A", "S", "M", "F", "G", "N"), col=colVec, las=1, ylim=c(2,yUpper))
abline(v=6.5, lty=2, lwd=.25, col="grey")
text(x=c(1:8), labels=rank(med), y=yUpper, col=colVec, cex=1.5)
dev.off()

#Plot testing losses
png(paste0(result_path, surv_model, "/testLoss_plot.png"))
med = apply(testLoss_mat, 2, median)
colfunc = colorRampPalette(c("chartreuse4","chartreuse2","yellow","orange","red"))
colVec = colfunc(8)[rank(med)]
par(cex.axis=1.25, cex.main=1.5) 
boxplot(as.list(as.data.frame(testLoss_mat)), main=paste0(toupper(surv_model), ": Testing Loss"), ylab="",
        names=c("B", "L", "A", "S", "M", "F", "G", "N"), col=colVec, las=1, ylim=c(2,yUpper))
abline(v=6.5, lty=2, lwd=.25, col="grey")
text(x=c(1:8), labels=rank(med), y=yUpper, col=colVec, cex=1.5)
dev.off()


############################### CREATE C-INDEX BOXPLOTS ##############################
#Read in training & testing c-indexes
trainConc_mat = read.table(paste0(result_path, surv_model, "/train_cIndex.csv"), sep=",", header=T)
testConc_mat = read.table(paste0(result_path, surv_model, "/test_cIndex.csv"), sep=",", header=T)

#Plot training c-indexes
png(paste0(result_path, surv_model, "/trainCIndex_plot.png"))
med = apply(trainConc_mat, 2, median)
colfunc = colorRampPalette(c("red","orange","yellow","chartreuse2","chartreuse4"))
colVec = colfunc(8)[rank(med)]
par(cex.axis=1.25, cex.main=1.5) 
boxplot(as.list(as.data.frame(trainConc_mat)), main=paste0(toupper(surv_model), ": Training C-Index"), ylab="",
        names=c("B", "L", "A", "S", "M", "F", "G", "N"), col=colVec, las=1, ylim=c(0.475,0.85))
abline(v=6.5, lty=2, lwd=.25, col="grey")
text(x=c(1:8), labels=rank(-med), y=0.85, col=colVec, cex=1.5)
dev.off()

#Plot testing c-indexes
png(paste0(result_path, surv_model, "/testCIndex_plot.png"))
med = apply(testConc_mat, 2, median)
colfunc = colorRampPalette(c("red","orange","yellow","chartreuse2","chartreuse4"))
colVec = colfunc(8)[rank(med)]
par(cex.axis=1.25, cex.main=1.5) 
boxplot(as.list(as.data.frame(testConc_mat)), main=paste0(toupper(surv_model), ": Testing C-Index"), ylab="",
        names=c("B", "L", "A", "S", "M", "F", "G", "N"), col=colVec, las=1, ylim=c(0.45,0.85))
abline(v=6.5, lty=2, lwd=.25, col="grey")
text(x=c(1:8), labels=rank(-med), y=0.85, col=colVec, cex=1.5)
dev.off()

