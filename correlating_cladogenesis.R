###DESCRIPTION###
#Based on these inputs, the code simulates the number and duration of sea level rises (or other periodic happenings that you are 
#interested in) that you've listed for the 'periods' variable (see `modify variables in this block` below) & randomly places these on 
#the time-scale back to your maxtimeline. It then counts how many of these simulated sea level rises coincide exactly with a
#cladogenesis event (the clade_ages variable for your tree) and writes this to "recordexact.txt". It also counts how many clades are 
#within 1mya of a simulated sea level rise, and writes this to recordrelaxed.txt. It redoes this simulation the number of times you list 
#in permutations.

#You can then pull these files into your favorite spreadsheet program and count how many of the simulations exceed your obsexact and 
#obsrelaxed values (it would also be really easy to add this into the script, but I like making a bar-graph of the distribution). This 
#is the p-value for for the number of cladogenesis events occurring at the same time as a sea level maxima, if they were in fact, not 
#correlated.

###MODIFY VARIABLES IN THIS BLOCK###
# how far back in time you want to go - this should be set to the earliest sea level inundation you have data for (that is older than your tree MRCA)
maxtimeline <- 30
# how fine scale you want the time divvied up - we were working in millions of years in this example
units <- 0.5
# The estimated MRCA of each clade in the tree
clade_ages <- c(8,4.5,15,3.5,9,9.5,11.5,17,24.5,16.5,4,2,2.5,5.5,6,3,3.5,4.5,5,6,12.5,18,19.5,20.5,28.5)
# The duration of each sea level fluctuation (or other periodic happening) listed separately for each event
periods <- c(0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,3,0.5)
# The number of clade MRCAs in the observed data found within a sea level maxima (or within periodic thing you are interested in).
obsexact <- 13
# The number of clade MRCAs in the observed data found within 1 mya of a sea level maxima (or within 1 mya of periodic thing you are interested in)
# This code doesn't actually use obsexact and obsrelaxed, but it would be easy to calculate the p-value using them.
obsrelaxed <- 21
# The number of simulations of your periodic things and how they correlate with your clade_ages that you want to perform
permutations <- 10000
###^^^MODIFY VARIABLES IN THIS BLOCK^^^###

timeline <- matrix("",ncol=2,nrow=((maxtimeline/units)-1))
timeline[,1] <- seq(1,maxtimeline,units)

timeline[(as.numeric(timeline[,1]) %in% clade_ages),2] <- "CLADE"

recordexact <- matrix(NA,ncol=1,nrow=permutations)
recordrelaxed <- matrix(NA,ncol=1,nrow=permutations)

for (j in 1:permutations) {

periodsample <- seq(1,maxtimeline,units)
temptimeline <- cbind(timeline,"")
temptimeline <- cbind(temptimeline,"")
checkout <- "no"

#Code review 15-Dec-2016: Not sure what the m=1000 if function will ever do, because I see no wway in the code I've incremented m to increase. I think this was a variable I included to troubleshoot the code when I first wrote it.
for (i in 1:(length(periods))) {
needaplace <- "yes"
m <- 1
while (needaplace=="yes") {
if (m==1000) {
checkout <- "yes"
break
}

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

if (checkout=="yes") {
break
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
