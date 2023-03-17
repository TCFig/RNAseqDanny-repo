# Align SRA reads to the Saccharomyces Cerevisiae genome
# copyright (c) 2022 - Danny Arends

#Create function "execute" - prevents downloading read files that you already have installed

execute <- function(x, outputfile = NA, intern = FALSE, quitOnError = FALSE){
  if(!is.na(outputfile) && file.exists(outputfile)){
    cat("Output for step exists, skipping this step\n");
    invisible("")
  }
  cat("----", x, "\n"); res <- system(x, intern = intern); cat(">>>>", res[1], "\n")
  if(res[1] >= 1){ 
    cat("Error external process did not finish\n\n");
    if(quitOnError) q("no")
  }
}

#Static variables
input.dir <- "~/data/raw"
input.base <- "SRR13978643" 
output.dir <- paste0("~/data/output/", input.base,".aln")
genome.path <- "~/genome/STAR/"
ref.fa.gz <- "~/genome/Saccharomyces_cerevisiae.R64-1-1.dna.primary_assembly.fa.gz"
ref.snps <- "~/genome/saccharomyces_cerevisiae.vcf.gz"

# Create an input folder
if(!file.exists(input.dir)){ dir.create(input.dir, recursive = TRUE) }

# Create an output folder
if(!file.exists(output.dir)){ dir.create(output.dir, recursive = TRUE) }


# STEP 0 - SRA Download and Compress
setwd(input.dir)

execute(paste0("fasterq-dump -p --split-files ", input.base), paste0(input.base, "_1.fastq"))
execute(paste0("bgzip ", input.base, "_1.fastq"), paste0(input.base, "_1.fastq.gz"))
execute(paste0("bgzip ", input.base, "_2.fastq"), paste0(input.base, "_2.fastq.gz"))

# STEP 1 - READ Trimming
trim.files  <- c(
                  paste0(input.dir, "/", input.base,"_1.fastq.gz"),
                  paste0(input.dir, "/", input.base,"_2.fastq.gz"),
                  paste0(output.dir, "/", input.base,"_1.P.fastq.gz"),
                  paste0(output.dir, "/", input.base,"_1.U.fastq.gz"),
                  paste0(output.dir, "/", input.base,"_2.P.fastq.gz"),
                  paste0(output.dir, "/", input.base,"_2.U.fastq.gz")
                )
trim.path <- "/home/tifigo/software/Trimmomatic"
trim.exec <- paste0("java -jar ", trim.path, "/dist/jar/trimmomatic-0.40-rc1.jar")
trim.opts <- paste0(" ILLUMINACLIP:", trim.path, "/adapters/TruSeq3-PE-2.fa:2:30:10")
trim.opts <- paste0(trim.opts, " LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36")
trim.cmd  <- paste0(trim.exec, " PE ", paste0(trim.files, collapse=" "), trim.opts)

execute(trim.cmd, trim.files[3])

# STEP 1.1 - UNZIP for STAR
execute(paste0("gunzip -k ", trim.files[3]), gsub(".fastq.gz", ".fastq", trim.files[3]))
execute(paste0("gunzip -k ", trim.files[5]), gsub(".fastq.gz", ".fastq", trim.files[5]))

files.in <- gsub(".fastq.gz", ".fastq", trim.files[c(3,5)])


# STEP 2 - Alignment using STAR
star.outbase <- paste0(output.dir, "/", input.base)
star.bam <- paste0(star.outbase, "Aligned.sortedByCoord.out.bam")

star.exec <- "STAR --runMode alignReads"
star.opts <- paste0("--genomeDir ", genome.path, " --outSAMtype BAM SortedByCoordinate")
star.in <- paste0("--readFilesIn ", paste0(files.in, collapse=" "))
star.out <- paste0("--outFileNamePrefix ", star.outbase)
star.cmd  <- paste0(star.exec, " ", star.in, " ", star.opts, " ", star.out)

execute(star.cmd, star.bam)

# STEP 2.1 - Create a samtools index
execute(paste0("samtools index ", star.bam), paste0(star.bam, ".bai"))

# STEP 2.2 - Create mapping and coverage statistics
execute(paste0("samtools flagstats ", star.bam))
execute(paste0("samtools coverage ", star.bam))


#STEP 3 - Remove duplicate reads using picard tools
p.bam <- paste0(star.outbase, "Aligned.sortedByCoord.RD.out.bam")
metrics.out <- paste0(star.outbase, "_metrics.txt")

p.exec <- "java -Xmx4g -jar ~/software/picard/build/libs/picard.jar"
p.in <- paste0("-I ", star.bam)
p.out <- paste0("-O ", p.bam, " -M ", metrics.out)
p.opts <- paste0("--REMOVE_DUPLICATES true")
p.cmd <- paste0(p.exec, " MarkDuplicates ", p.opts," ", p.in, " ", p.out)

execute(p.cmd, p.bam)

# STEP 3.1 - Create a samtools index
execute(paste0("samtools index ", p.bam), paste0(p.bam, ".bai"))

# STEP 3.2 - Create mapping and coverage statistics
execute(paste0("samtools flagstats ", p.bam))
execute(paste0("samtools coverage ", p.bam))


# STEP 4 - Add read group (1) and sample run, library, and name
rg.bam <- paste0(star.outbase, "Aligned.sortedByCoord.RD.RG.out.bam")
rg.opts <- paste0("-PL ILLUMINA -PU run -LB ", gsub("SRR", "", input.base), " -SM ", input.base)
p.cmd <- paste0(p.exec, " AddOrReplaceReadGroups -I ", p.bam, " -O ", rg.bam, " ", rg.opts)
execute(p.cmd)

# STEP 4.1 - Create a samtools index
execute(paste0("samtools index ", rg.bam), paste0(rg.bam, ".bai"))


# STEP 5 - GATK prep
gatk.exec <- "java -Xmx4g -jar ~/software/gatk-4.3.0.0/gatk-package-4.3.0.0-local.jar"
gatk.opts <- paste0("-R ", ref.fa.gz, " --known-sites ", ref.snps)

# STEP 5.1 - GATK BaseRecalibrator
gatk.cov1 <- paste0(star.outbase, "_cov1.txt")
gatk.cmd  <- paste0(gatk.exec, " BaseRecalibrator ", gatk.opts, " -I ", rg.bam, " -O ", gatk.cov1)
execute(gatk.cmd, gatk.cov1)

# STEP 5.2 - GATK ApplyBQSR
recal.bam <- paste0(star.outbase, "Aligned.sortedByCoord.RD.RG.RC.out.bam")
gatk.cmd  <- paste0(gatk.exec, " ApplyBQSR -R ", ref.fa.gz, " -bqsr ", gatk.cov1, " -I ", rg.bam, " -O ", recal.bam)
execute(gatk.cmd, recal.bam)

# STEP 5.3 - GATK BaseRecalibrator
gatk.cov2 <- paste0(star.outbase, "_cov2.txt")
gatk.cmd  <- paste0(gatk.exec, " BaseRecalibrator ", gatk.opts, " -I ", recal.bam, " -O ", gatk.cov2)
execute(gatk.cmd, gatk.cov2)

# STEP 5.4 - GATK AnalyzeCovariates
recal.plot <- paste0(star.outbase, "AnalyzeCovariates.pdf")
gatk.cmd  <- paste0(gatk.exec, " AnalyzeCovariates -before ", gatk.cov1, " -after ", gatk.cov2, "  -plots ", recal.plot)
execute(gatk.cmd)


# STEP 6 - Index the recalibrated bam files
execute(paste0("samtools index ", recal.bam), paste0(recal.bam, ".bai"))

# STEP 6.1 - Create mapping and coverage statistics
execute(paste0("samtools flagstats ", recal.bam))
execute(paste0("samtools coverage ", recal.bam))

q("no")
