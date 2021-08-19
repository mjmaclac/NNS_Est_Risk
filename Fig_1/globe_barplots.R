# Import data, save out a reshaped wide format for use in Arc, long format works fine for R and QGIS?
rm(list = ls())                                               # Clears previous information
dev.off()                                                     # Closes plots    
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))   # Sets working directory to this file's location

source("Preamble.R")  

# Imports and species discovery data
imports <- read.csv("Data/Imports-data-v4.csv", stringsAsFactors = FALSE)
disc    <- read.csv("Data/Discoveries-data-v7.csv", stringsAsFactors = FALSE)

# Remove total and Neararctic regions
imports <- imports[!(imports$Region %in% c('Total', 'Nearctic')), ]
disc    <- disc[!(disc$Region %in% c('Total', 'ROW', 'Nearctic','Aggregate')), ]

# Drop years in discovery data before 1854 to match trade data
disc <- disc[!(disc$Year < 1800), ]

# Get species totals by region for all years.
polylabsAll <- disc %>% group_by(Region) %>% summarise_each(funs(sum))  

disc <- disc[!(disc$Year < 1854), ]

# Factor variable to bin values (every 5 years, except 1854-1857 (4 years))
imports$indx <- factor(c(rep(1,4), rep(2:32, each=5))) 
disc$indx    <- factor(c(rep(1,4), rep(2:32, each=5)))

impbin <- imports %>% group_by(Region, indx) %>% summarize(Value = sum(Value)) 

discbin <- disc %>% group_by(Region, indx) %>% summarize(Discoveries = sum(Discoveries)) 

# Start of period (bin)
x    <- 1857
m    <- rep(NA, 32)
n    <- rep(NA, 32)
m[1] <- 1854
n[1] <- 1857
for (i in 1:31) {
  x      <- x + 5
  n[i+1] <- x
  m[i+1] <- n[i+1] - 4
}

levs       <- data.frame(start = m, end = n)
levs$label <- paste(as.character(levs$start), as.character(levs$end), sep = '-')

# Region buffers
rb <- shapefile("Regional Inset Buffer (150km) Shapefiles/Merged Region Inset Buffer Shapefile/Region_Inset_Buffers.shp")

# Points
p <- shapefile("Region Points Shapefile/Region_Points.shp") 
# filter
p <- p[!(p$Region=="Nearctic"),]
# sort
p <- p[order(p$Region),]

# Flow lines
fl <- shapefile("Flow Lines Shapefile/Region_FlowLines.shp")
fln <- fl[order(fl$Region),]

# Global and US shapefiles
shp <- shapefile("Buffer Shapefiles/Clippings/Global Regions Shapefiles/Global_Regions_Dissolved.shp")
us_shp <- shapefile("United States Shapefile/US_withAlaska.shp")

# For plotting
cols    <- rev(c("#f57ab6", "#018a71", "#fa8500", "#87b833", "#fa4632", "#0073ff", "#956cb2"))
theme1  <- c("#018a71", "#fa4632", "#fa8500", "#956cb2", "#87b833", "#cfba00", "#0073ff", "#f57ab6")
flncols <- c("#018a71", "#fa4632", "#fa8500", "#956cb2", "#87b833", "#0073ff", "#f57ab6")

regions <- unique(imports$Region)

# Wide-format matrix for barplot
cimports <- spread(impbin, indx, Value)
cimports <- as.matrix(cimports[-1])
rownames(cimports) <- regions
colnames(cimports) <- levs$label

# Wide-format matrix for barplot
cdisc <- spread(discbin, indx, Discoveries)
cdisc <- as.matrix(cdisc[-1])
rownames(cdisc) <- regions
colnames(cdisc) <- levs$label


# Want to sort by cummulative imports
ord <- names(rev(sort(rowSums(cimports))))
cdisc <- cdisc[ord, ]
cimports <- cimports[ord, ]

# For a gray background
bluemarble_US <- brick("BlueMarble2004-naclean_USRobison.tif")
r <- bluemarble_US[[1]]
r[!is.na(r)] <- 1 
er <- extent(-16986823, 16992938, -6500000, 8615716)
rc <- crop(r, er)
rca <- aggregate(rc, fact = 3)

pol <- rasterToPolygons(rca, dissolve = TRUE)  #takes time...30 sec?


# For polygons labels (cummulative sum of species discovery )
polylabs <- rowSums(cdisc)[c(6, 3, 5, 1, 4, 2, 7)]
polylabs[1:7] <- polylabsAll$Discoveries  # replace with cumulative sum for all years.
names(polylabs)[4] <- "European Palearactic"
pl <- data.frame(Region = names(polylabs), cumsum = polylabs)
shp <- merge(shp, pl, by = 'Region')



windows(14,10)
plot(shp)
plabs <- polygonsLabel(shp[!is.na(shp$Region) & shp$Region != 'Nearctic', ], 
                       labels = shp$cumsum[!is.na(shp$cumsum)], 
                       method = 'centroid', col = 'white', cex = 3, doPlot = FALSE)
plabs[1,1] <- 11450000
plabs[2,1] <- -14500000
plabs[4,2] <- 5700000
plabs[6,1] <- 3501138
plabs[5,1] <- -12300000
plabs[5,2] <- 500000
plabs[7,1] <- -8500000
plabs[7,2] <- -2000000

# For lines showing cummulative import value
rsums <- rowSums(cimports)                   # Get cumulative imports for each region for entire history
linewidth <- (100 * (log(rsums + 1))) / 8    # Scale those values to make for good lines for figure.
names(linewidth)[1] <- "European Palearactic"  

lw <- data.frame(Region = names(linewidth), width = linewidth)
fln <- merge(fln, lw, by = 'Region')

# Plot
png("./Figure_1.png",width = 10000, height = 8000, units = "px", res = 800)
  layout(matrix(c(rep(1,6),rep(2,4)), 10, 1, byrow = TRUE))
  
  par(mar = c(8,11,2,14))
  
  # Map
  plot(pol, col = 'lightblue2', border = 'transparent')
  plot(shp[!is.na(shp$Region), ], col=theme1, add = TRUE)
  plot(fln[!(fln$Region == "Nearctic"),], lwd = fln$width, col = flncols, add = TRUE)
  plot(shp[!is.na(shp$Region), ], col = theme1, add = TRUE)
  text(plabs[,1], plabs[,2], labels = shp$cumsum[!is.na(shp$cumsum)], cex = 3, col = 'white')
  plot(us_shp, col = "#ebd628", add = TRUE)
  
  #create barplot for imports
  par(new=TRUE)
  par(mar=c(4,10,3,1), mgp = c(6,2,0))
  
  bp <- barplot(cimports, col = cols, axes = FALSE, axisnames=FALSE, border = NA,
                ylab = "Nursery products import value ($2015B)", cex.lab = 2.25, 
                ylim = c(0, 3))
  axis(2, cex.axis = 2, las = 1, tck = -.02)
  
  mtext("(A)", side=2, line=6, padj=-10, cex=2, las=2, font=2)  # side 2 is left, padj sets vertical, las 2 sets horizonatl orientation, font = 2 sets bold

  par(mar = c(8,10,0,1), mgp = c(6,2,0))
  
  labs <- c("1854-1857", "1888-1892", "1928-1932", "1968-1972", "2008-2012")  # 
  
  bp2 <- barplot(cdisc, col = cols, axes = FALSE, axisnames = FALSE, border = NA, 
                 ylab="Species discovered", xlab = "Years", cex.lab = 2.5, ylim = c(0, 50))
  legend(-1,52, legend = shp$Region[!is.na(shp$Region)], col=theme1, border="black", bty="o", 
         bg=NULL, box.lty=0, cex=2, fill=theme1, ncol=1, y.intersp=1)
  mtext("(B)", side=2, line=6, padj=-6, cex=2, las=2, font=2)  # side 2 is left, padj sets vertical, las 2 is horizontal orientation
  staxlab(1, at = c(bp[1], bp[8], bp[16], bp[24], bp[32]), labels = labs, cex = 1.5, top.line = 1.78, line.spacing = 0.01)
  axis(2, cex.axis = 2, las = 1, tck = -.04)  
dev.off()
