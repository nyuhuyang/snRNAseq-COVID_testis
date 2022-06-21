Samples = list.files("data/fastq")
setwd("data/fastq")
rm.names = paste(c("BRI1469A13_","BRI1469A19_","BRI1469A15_",
                   "BRI1469A11_","BRI1469A16_","BRI1469A14_",
                   "BRI1469A20_","BRI1469A18_","BRI1469A12_",
                   "BRI1469A17_"),collapse = "|")
for(i in seq_along(Samples)){
    fastq.files = list.files(Samples[i])
    for(file in fastq.files){
        file.rename(from = paste0(Samples[i],"/",file),
                    to = paste0(Samples[i],"/",sub(rm.names,"",file)))
    }
    print(list.files(Samples[i]))
}