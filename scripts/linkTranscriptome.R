library(tximeta)
indexDir <- file.path("data/airway", "gencode.v43.salmon")
fasta <- "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_43/gencode.v43.transcripts.fa.gz"
gtf <- "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_43/gencode.v43.chr_patch_hapl_scaff.basic.annotation.gtf.gz"
makeLinkedTxome(indexDir=indexDir,
                source="LocalGencode",
                organism="Homo sapiens",
                release="43",
                genome="chr38",
                fasta=fasta,
                gtf=gtf,
                write=FALSE)



