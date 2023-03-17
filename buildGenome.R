# Download Saccharomyces Cerevisiae genome
# copyright (c) 2022 - Danny Arends
#

uri <- "ftp.ensembl.org/pub/release-109/fasta/saccharomyces_cerevisiae/dna/"
base <- "Saccharomyces_cerevisiae.R64-1-1.dna.chromosome."
chrs <- c(as.character(as.roman(seq(1:16))), "Mito")

# Download
for(chr in chrs){
  fname <- paste0(base, chr, ".fa.gz")

  # Download command
  cmd <- paste0("wget ", uri, fname)
  #cat(cmd, "\n")
  system(cmd)
}


# Create an empty the file
cat("", file = "Saccharomyces_cerevisiae.R64-1-1.dna.primary_assembly.fa")
for(chr in chrs){
  fname <- paste0(base, chr, ".fa.gz")
  # Extract and merge into a fast file
  cmd <- paste0("zcat ", fname, " >> Saccharomyces_cerevisiae.R64-1-1.dna.primary_assembly.fa")
  #cat(cmd, "\n")
  system(cmd)
}

# Compress the fasta file using bgzip (keep original)
cmd <- paste0("bgzip -k Saccharomyces_cerevisiae.R64-1-1.dna.primary_assembly.fa")
#cat(cmd, "\n")
system(cmd)

# Delete the chromosomes
for(chr in chrs){
  fname <- paste0(base, chr, ".fa.gz")
  # Extract and merge into a fast file
  cmd <- paste0("rm ", fname)
  #cat(cmd, "\n")
  system(cmd)
}
