############
# Preamble #
############

# Packages to be loaded 
packages <- c("dplyr", "tidyr","RODBC", "xlsx","blsAPI","lubridate",
              "zoo","forecast","tseries","stringr","ggplot2","egg","ggpubr",
              "rgeos","data.table","RColorBrewer","raster","plotrix")

# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

# Packages loading
invisible(lapply(packages, library, character.only = TRUE))