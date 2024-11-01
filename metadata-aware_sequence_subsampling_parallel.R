library(ape)
library(Biostrings)
library(parallel)
library(foreach)

# Set up parallelization
n.cores <- parallel::detectCores() - 1
# Create the cluster
my.cluster <- parallel::makeCluster(n.cores, type = "PSOCK")
print(my.cluster)
suppressWarnings(doParallel::registerDoParallel(cl = my.cluster))


l1 <- read.csv("l1_sequences.tsv", header = T, sep = "\t")
#Check country and host info just to make sure things make sense
table(l1$Country)
table(l1$Host)
#Filter empty entries in Host and Country columns
l1 <- l1[!(is.na(l1$Host) | l1$Host==""),]
l1 <- l1[!(is.na(l1$Country) | l1$Country==""),]
#Set up vectors of countries to define geo regions
North_America <- c("Canada", "USA")
South_America <- c("Brazil", "Chile")
Africa <- c("Cameroon", "South Africa")
Asia <- c("China", "Indonesia", "Malaysia", "South Korea")
Europe <- c("Germany", "Hungary", "Norway", "Slovenia")
#Go through the entries in the Country column and replace country with a corresponding region
for(i in 1:length(l1$Country)){
    if(l1$Country[i] %in% North_America) {l1$Country[i] = "North_America"}
    if(l1$Country[i] %in% South_America) {l1$Country[i] = "South_America"}
    if(l1$Country[i] %in% Asia) {l1$Country[i] = "Asia"}
    if(l1$Country[i] %in% Africa) {l1$Country[i] = "Africa"}
    if(l1$Country[i] %in% Europe) {l1$Country[i] = "Europe"}
}
#Check that lumping of countries into geo regions went fine and there is no ungrouped country left
table(l1$Country)

#Get all accession numbers in the dataset
accessions <- as.vector(l1$Accession)

#Pull sequence records from GenBank
x <- read.GenBank(accessions[1:length(accessions)], as.character = T)

toDiscard <- c()

toDiscard <- foreach(i=1:length(x), .combine='c', .init=c(), .multicombine = T) %dopar% {
    print(i)
    s1 = gsub(", ", "", toString(x[[i]]))
    for(j in i+1:length(x){
        if(j <= length(x)){
            s2 = gsub(", ", "", toString(x[[j]]))
            palign1 <- Biostrings::pairwiseAlignment(s1, s2)
            identity = Biostrings::pid(palign1)
            if(identity > 99.5 && l1$Country[i] == l1$Country[j] && l1$Host[i] == l1$Host[j]){
                toDiscard <- c(toDiscard, l1$Accession[i])
            }
        }
    }
return(unique(toDiscard))
}
stopCluster(my.cluster)
#Assemble a list of accession numbers to keep (deduplicated list) and save to file
keep <- setdiff(accessions, toDiscard)
keep_seqs <- read.GenBank(keep, as.character = T)
write.dna(keep_seqs, "l1_deduplicated.fasta", format = "fasta")