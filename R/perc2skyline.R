# Convert Percolator txt output to ssl file for Skyline MS1 filtering
library(tidyverse)

percDat <- read_tsv("data/crux-output/percolator.target.psms.txt")

fileKey <- tibble(file = c("180405-BTAP-LRP-Rep1-0.mzML",
                           "180405-BTAP-LRP-Rep1-1.mzML",
                           "180405-BTAP-LRP-Rep1-2.mzML",
                           "180405-BTAP-LRP-Rep1-3.mzML",
                           "180405-BTAP-LRP-Rep1-4.mzML",
                           "180405-BTAP-LRP-Rep1-5.mzML",
                           "180405-BTAP-LRP-Rep1-6.mzML",
                           "180405-BTAP-LRP-Rep1-7.mzML",
                           "180405-BTAP-LRP-Rep1-GST.mzML",
                           "180405-BTAP-LRP1B-Rep1-0_180427065706.mzML",
                           "180405-BTAP-LRP1B-Rep1-1.mzML",
                           "180405-BTAP-LRP1B-Rep1-2.mzML",
                           "180405-BTAP-LRP1B-Rep1-3.mzML",
                           "180405-BTAP-LRP1B-Rep1-4.mzML",
                           "180405-BTAP-LRP1B-Rep1-5.mzML",
                           "180405-BTAP-LRP1B-Rep1-6.mzML",
                           "180405-BTAP-LRP1B-Rep1-7.mzML",
                           "180417-BTAP-LRP1B-Rep1-GST-Scouting.mzML"),
                  file_idx = c(0:(length(file)-1)))

sslOut <- percDat %>%
    left_join(fileKey) %>%
    mutate(`score-type` = "PERCOLATOR QVALUE") %>%
    select(file, scan, charge, sequence,`score-type`, `percolator q-value`) %>%
    rename(score = `percolator q-value`) %>%
    mutate(sequence = str_replace(sequence, "15.99", "+15.99")) %>%
    filter(score <= 0.01)

write_tsv(sslOut, "data/percolator_results.ssl")
