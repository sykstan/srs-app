
require(data.table)
require(magrittr)
require(tidyr)

sbp.ccsdt <- read.csv(file = "~/GoogleDrive/srs-app/data/sbp-ccsdt-adz.edat", header = TRUE
         , sep = "|", strip.white = TRUE) %>% data.table()