## Started by Alice Wang to show ESR1 wave.
##

require(bigWig)
step.bpQuery.bigWig
load.bigWig

windowsize <- 3000


B70 <- load.bigWig("../data/MCF-7_B7_GDNF_0_hr_plus.bw")
ESRB70 <- step.bpQuery.bigWig(B70,"chr6", 152128453, 152424409, windowsize)
B70mean <- abs(B70$mean)
B70bc <- B70$basesCovered
B70reads <- B70mean*B70bc

B71 <- load.bigWig("../data/MCF-7_B7_GDNF_1_hr_plus.bw")
ESRB71 <- step.bpQuery.bigWig(B71,"chr6", 152128453, 152424409, windowsize)
B71mean <- abs(B71$mean)
B71bc <- B71$basesCovered
B71reads <- B71mean*B71bc


Difference <- ESRB71/B71reads-ESRB70/B70reads
Bases <- c(1:NROW(Difference))*windowsize-(windowsize/2)
plot(Bases, Difference, cex=2, pch=19)
abline(h=0,col='red')

#Calculating end of wave 
require(groHMM)

source("polymeraseWave.bw.R")
ESR <- data.frame(chr="chr6", start=152128453, end=152424409, str="+", genes="ESR1", score=0)
elongationRate <- polymeraseWaveBW(
"/local/storage/projects/mcf7tamres/data/MCF-7_B7_GDNF_0_hr_plus.bw",
"/local/storage/projects/mcf7tamres/data/MCF-7_B7_GDNF_0_hr_minus.bw",
"/local/storage/projects/mcf7tamres/data/MCF-7_B7_GDNF_1_hr_plus.bw",
"/local/storage/projects/mcf7tamres/data/MCF-7_B7_GDNF_1_hr_minus.bw",
ESR, approxDist=50000, returnVal="alldata", prefix="IMG/ESR1.", emissionDistAssumption= "gamma", size=1000)

60- (104 / 2.1) ## Start time, assuming median elongation rate in MCF-7 cells.

## Elongation rate in MCF-7 cells.
wave40 <- polymeraseWaveBW(
"/local/storage/data/hg19/mcf7/groseq/MCF7.unt.all_plus.bw",
"/local/storage/data/hg19/mcf7/groseq/MCF7.unt.all_minus.bw",
"/local/storage/data/hg19/mcf7/groseq/MCF7.40m.all_plus.bw",
"/local/storage/data/hg19/mcf7/groseq/MCF7.40m.all_minus.bw",
ESR, approxDist=20000, returnVal="alldata", prefix="IMG/ESR1-Hah40.", emissionDistAssumption= "gamma", size=1000)

wave10 <- polymeraseWaveBW(
"/local/storage/data/hg19/mcf7/groseq/MCF7.unt.all_plus.bw",
"/local/storage/data/hg19/mcf7/groseq/MCF7.unt.all_minus.bw",
"/local/storage/data/hg19/mcf7/groseq/MCF7.10m.all_plus.bw",
"/local/storage/data/hg19/mcf7/groseq/MCF7.10m.all_minus.bw",
ESR, approxDist=50000, returnVal="alldata", prefix="IMG/ESR1-Hah10.", emissionDistAssumption= "gamma", size=1000)


## Estimate start time: 
60- (104/ ((wave40[[2]]$EndWave-wave10[[2]]$EndWave)/1000/30))


