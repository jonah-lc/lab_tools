

#### In this example, 2 titrations were performed to assess the binding of the DNA reported to POPC
#### LUVs extruded to two 200 and 50 nm. Each titration was ran in duplicate on an agarose gel

#### Band intensities were extracted using ImageJ, the data is included. 
#### Each csv is 14 rows, the first 7 rows are the band intensities for the titration. Lane 1 
#### is the 0% bound control (ie no LUVs) and lane 7 is the 100% bound control (ie no DNA)
#### The next 7 lanes are the background intensities for each of the first 7 lanes


# import tools
source("/media/jonah/Lickitung/SYNC/4-SCRIPTS/R/jlc_densitomtery_tools.R")


## set labels for plot
xlab = "\n[Vesicles] / nM"
ylab = "% Bound\n"
labs = c("POPC-200", "POPC-50")
## set colours for plot
cols = c("#ff8253","#7d270b")

## set labels for data
## data directory, you will need to adjust this
dir = "/media/jonah/Lickitung/SYNC/8-LAB_DATA/TITRATION/20211119_POPC&Fuse/"
## file names
popc2 = c("popc200-1", "popc200-2")
popc5 = c("popc50-1", "popc50-2")
## extension types
ext =".csv"

## concentration of LUVs used (uM)
popc_uM = c(0, 25, 62.5, 125, 250, 875)
## convert to concentration of vesicles
## hydrodynamic diamateres were assessed by DLS
## conc_mem2particle needs concentrations in nM. The default parameters are for POPC vesicles
## for other vesicles you may need to adjust the head group area
nM_p200 = conc_mem2particle(popc_uM/1000, hydro_diam = 180)
nM_p50 = conc_mem2particle(popc_uM/1000, hydro_diam = 70)

## nM vesicle values to fit model between
startval = 0.001
endval = 200

## prepare data 
p200  = densCalc(dir, popc2, ext, nM_p200, labs[1])
p50 = densCalc(dir, popc5, ext, nM_p50, labs[2])

## build models
mp200 = fitter(p200, labs[1], 1000, startval, endval)
mp50 = fitter(p50, labs[2], 1000, startval, endval)

## stack data frames for easier plotting
df = rbind(p200, p50)
mf = rbind(mp200, mp50)

## plot
plt = densPlot(df, mf, cols, xlab, ylab, logax="x")
plt

## save plot
pngOut(plt, "output/diretory/of/your/choice")

## extract Kd
print(paste0(unique(mf$n_id),", Kd = ", unique(mf$km), " +/- ",unique(mf$error)))
