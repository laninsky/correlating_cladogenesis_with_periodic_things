# correlating_cladogenesis_with_periodic_things

````
maxtimeline <- 30
units <- 0.5
clade_ages <- c(8,4.5,15,3.5,9,9.5,11.5,17,24.5,16.5,4,2,2.5,5.5,6,3,3.5,4.5,5,6,12.5,18,19.5,20.5,28.5)
periods <- c(0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,3,0.5)
obsexact <- 13
obsrelaxed <- 21
permutations <- 10000

timeline <- matrix("",ncol=2,nrow=((maxtimeline/units)-1))
timeline[,1] <- seq(1,maxtimeline,units)

timeline[(as.numeric(timeline[,1]) %in% clade_ages),2] <- "CLADE"

recordexact <- matrix(NA,ncol=1,nrow=permutations)
recordrelaxed <- matrix(NA,ncol=1,nrow=permutations)

for (j in 1:permutations) {

periodsample <- seq(1,maxtimeline,units)
temptimeline <- cbind(timeline,"")
temptimeline <- cbind(temptimeline,"")

for (i in 1:(length(periods))) {
needaplace <- "yes"
while (needaplace=="yes") {
x <- sample(periodsample, 1, replace = TRUE, prob = NULL)

if (!((periods[i]+x)>maxtimeline)) {
no_of_slots <- periods[i]/units
if(no_of_slots==1) {
if (temptimeline[(as.numeric(timeline[,1])==x),3]=="") {
temptimeline[(as.numeric(timeline[,1])==x),3] <- "PERIOD"
needaplace <- "no"
periodsample <- periodsample[periodsample!=x] 
}
} else {
rowno <- which(as.numeric(timeline[,1])==x)
if(all(temptimeline[rowno:(rowno+no_of_slots),3]=="")) {
for (y in 0:no_of_slots) {
temptimeline[(rowno+y),3] <- "PERIOD"
}
needaplace <- "no"
periodsample <- periodsample[periodsample!=x] 
}
}
}
}
}

if (is.null(dim(temptimeline[(temptimeline[,2]=="CLADE" & temptimeline[,3]=="PERIOD"),])[1])) {
recordexact[j] <- 0 } else {
recordexact[j] <- dim(temptimeline[(temptimeline[,2]=="CLADE" & temptimeline[,3]=="PERIOD"),])[1]
}

for (k in 1:((maxtimeline/units)-1)) {
if(temptimeline[k,3]=="PERIOD") {
if (k>2 & (k < ((maxtimeline/units)-2))) {
temptimeline[k,4] <- "PERIOD"
temptimeline[(k-2),4] <- "PERIOD"
temptimeline[(k-1),4] <- "PERIOD"
temptimeline[(k+1),4] <- "PERIOD"
temptimeline[(k+2),4] <- "PERIOD"
} else {
if (k <= 2) {
temptimeline[1:k,4] <- "PERIOD"
temptimeline[(k+1),4] <- "PERIOD"
temptimeline[(k+2),4] <- "PERIOD"
} else {
temptimeline[(k-2),4] <- "PERIOD"
temptimeline[(k-1),4] <- "PERIOD"
temptimeline[k:((maxtimeline/units)-1),4] <- "PERIOD"
}
}
}
}

if (is.null(dim(temptimeline[(temptimeline[,2]=="CLADE" & temptimeline[,4]=="PERIOD"),])[1])) {
recordrelaxed[j] <- 0 } else {
recordrelaxed[j] <- dim(temptimeline[(temptimeline[,2]=="CLADE" & temptimeline[,4]=="PERIOD"),])[1]
}
}

write.table(recordrelaxed,"recordrelaxed.txt")
write.table(recordexact,"recordexact.txt")
