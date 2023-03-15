library(Biostrings)
library(ShortRead)


x <- readDNAStringSet("Documents/Agilent/customers/Interpace/202302/Project10/WES-230223-36408271-MEGL-MCH-7_SUBSET_S7_L001_R1_001.fastq.gz")


x <- readFastq("Documents/Agilent/customers/Interpace/202302/Project10/WES-230223-36408271-MEGL-MCH-7_SUBSET_S7_L001_R1_001.fastq.gz")

head(sread(x))
head(quality(x))

encoding(quality(x))


fq <- fs::dir_ls("Documents/Agilent/customers/Interpace/202302/Project9/", regexp = "fastq.gz$")
qa <- qa(fq, type="fastq")
browseURL(report(qa))

ll <- qa[["perCycle"]]$quality$lane == "nicolo-AT-RLC21-1-z82_S267_R2_001.fastq.gz"

ShortRead:::.plotCycleQuality(qa[["perCycle"]]$quality[ll,])


qq <- qa[["perCycle"]]$quality[ll,]

library(tidyverse)

qq |> 
  group_by(Cycle) |> 
  filter(Count != min(Count)) |> 
  mutate(Median = weighted.mean(Score, Count)) |> 
  ggplot(aes(x = Cycle,
             y = Median)) +
  geom_line() +
  ylim(0, 40)


qq |> 
  count(Cycle)


qq |> 
  group_by(Cycle) |>
  mutate(somma= sum(Count))