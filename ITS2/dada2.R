library(dada2)

# Should have all sequencing files stored in a folder called Cornell16S
# The dada2.sh script will need to unzip these before running R script
path = "./Cornell16S"
list.files(path)

# Sort ensures forward/reverse reads are in same order
fnFs = sort(list.files(path, pattern="_R1_001.fastq.gz"))
fnRs = sort(list.files(path, pattern="_R2_001.fastq.gz"))


# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names = sapply(strsplit(fnFs, "L001_"), `[`, 1)

# Specify the full path to the fns
fnFs = file.path(path, fnFs)
fnRs = file.path(path, fnRs)

# Place filtered files in filtered/ subdirectory
filt_path = file.path(path, "filtered")
filtFs = file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs = file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))


# Filter and trim sequence files - error rates are high 
out = filterAndTrim(fnFs, filtFs, fnRs, filtRs,  trimLeft=c(12,12),truncLen=c(200,170),maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE, n=1e6,verbose=TRUE ,compress=TRUE, multithread=TRUE)
# Used maxEE=2 and trucQ=2 - defaults
# Will want to try different cutoffs to see if that helps with quality issues

# Update the paths to files to include only those that didn't throw out all seqs
filtFs = file.path(filt_path, sort(list.files(filt_path, pattern="_F_filt.fastq.gz")))
filtRs = file.path(filt_path, sort(list.files(filt_path, pattern="_R_filt.fastq.gz")))

# Update the sample names to include only those that didn't throw out all seqs
derepFs.names = sapply(strsplit(sort(list.files(filt_path, pattern="_F_filt")), ".assembled"), `[`, 1)
derepRs.names = sapply(strsplit(sort(list.files(filt_path, pattern="_R_filt")), ".assembled"), `[`, 1)


# Learn errors - using all reads because data are sparse
errF = learnErrors(filtFs, nreads=4e+06, randomize=TRUE, multithread=TRUE)
errR = learnErrors(filtRs, nreads=4e+06, randomize=TRUE, multithread=TRUE)

# Making plots
pErrF = plotErrors(errF, nominalQ=TRUE)
pErrR = plotErrors(errR, nominalQ=TRUE)

# Dereplicating sequences to avoid duplicates
derepFs = derepFastq(filtFs, 4e+06)
derepRs = derepFastq(filtRs, 4e+06)

# Naming the dereplicated sequences
names(derepFs) = derepFs.names
names(derepRs) = derepRs.names

# Running dada2 on seqs
dadaFs = dada(derepFs, err=errF, selfConsist=FALSE, pool=TRUE,OMEGA_A = 1e-40, USE_QUALS=TRUE,BAND_SIZE=16,VERBOSE=FALSE, multithread=TRUE)
dadaRs = dada(derepRs, err=errR, selfConsist=FALSE, pool=TRUE,OMEGA_A = 1e-40, USE_QUALS=TRUE,BAND_SIZE=16,VERBOSE=FALSE, multithread=TRUE)
# First full run default settings OMEGA_A = 1e-40 and BAND_SIZE=16

mergers = mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
# Merging dereplicated sequences, default values require no mismatches
# minOverlap default is 20


# Creating OTU table
OTUtab = makeSequenceTable(mergers)

# Removing chimeras - reiterating after first dada step, based on merged seqs
OTUtab.nochim = removeBimeraDenovo(OTUtab, method="consensus", multithread=TRUE, verbose=TRUE)

# Getting the sequence losses along the way
getN = function(x) sum(getUniques(x))
track = cbind(out, sapply(mergers, getN), rowSums(OTUtab), rowSums(OTUtab.nochim))
colnames(track) = c("input", "filtered", "denoised", "tabled", "nonchim")
rownames(track) = sample.names
head(track)

# Saving the OTU sequences as a fasta file
uniquesToFasta(OTUtab.nochim,fout="DADA2_seqs_nochim.fasta")
uniquesToFasta(OTUtab,fout="DADA2_seqs.fasta")

# Saving the OTU table as R object in the local directory
saveRDS(OTUtab.nochim, "OTUtab.nochim.rds")
saveRDS(OTUtab, "OTUtab.rds")
saveRDS(track, "track.rds")
saveRDS(pErrF, "pErrF.rds")
saveRDS(pErrR, "pErrR.rds")

# Need to have the .sh script bring out OTUtab.rd, DADA2_seqs.fasta, track.rds, pErrR.rds and pErrF.rds -zip them all together, and save that to Gluster.
