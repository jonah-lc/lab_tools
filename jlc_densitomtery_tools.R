#### Functions for plotting densitometry of gels #####se



source("/media/jonah/Lickitung/SYNC/4-SCRIPTS/R/jlc_lab_tools.R")

library(drc)
library(scales)

densImport = function(path, id) {
  ## imports data set, splits to 2x columns, sample and background
  ## path   (str)     absolute path to data files
  ## id     (str)     identifier
  dat = read.csv(path)
  df_len = nrow(dat)
  df = data.frame(dat[1:(df_len/2),], dat[((df_len/2)+1):nrow(dat),])
  names(df) = c("ints", "bkg")
  df = na.omit(df)
  df$id = sub(".*/", "", path)
  return(df)
}


perc_bound = function(imported_df){
  # calculates percentage binding from dens data
  ## imported_df    (df)    output from "densImport" module
  
  ## background correction
  imported_df$corrected = imported_df$ints - imported_df$bkg
  ## normalise
  imported_df$norm = norm(imported_df$corrected)
  ## turn to percentage
  imported_df$pb = (1-imported_df$norm)*100
  return(imported_df)
}



preAvg = function(folder, file, ext=".csv", binder_conc){
  inp = densImport(paste0(folder, file, ext))
  inp = perc_bound(inp)
  ## now remove last line post normalisation
  inp = inp[1:(nrow(inp)-1), ]
  ## add new titratant concs
  inp$binder = binder_conc

return(inp)
}
  
comSub<-function(x) {
  # Not my function, taken from https://stackoverflow.com/questions/28273716/r-implementation-for-finding-the-longest-common-starting-substrings-in-a-set-of
  
  ## will find first common substring of any list of strings
  # sort the vector
  x<-sort(x)
  # split the first and last element by character
  d_x<-strsplit(x[c(1,length(x))],"")
  # search for the first not common element and so, get the last matching one
  der_com<-match(FALSE,do.call("==",d_x))-1
  # if there is no matching element, return an empty vector, else return the common part
  ifelse(der_com==0,return(character(0)),return(substr(x[1],1,der_com)))
}



densCalc= function(folder, file_list, ext=".csv", binder_conc, label){
  # binder_conc     (scalar)      list of int for concentration of titrator, must be same half length of data files
  
  ## prepare data
  df = do.call(rbind, lapply(file_list, preAvg, folder=folder, ext=ext, binder_conc))
  # how many repeats
  dat_numb = length(unique(df$id))
  ## avg and sd
  avg = aggregate(pb ~ binder, df, mean)
  stdev = aggregate(pb ~ binder, df, sd)
  # put into df
  out  = data.frame(x=binder_conc, y=avg[,2], stdev=stdev[,2], id=label)
  return(out)
}
  



fitter = function(df, id, length_out=1000, startval = 0.001, endval=100){
  ## will fit model to data from dens_cals.
  ## df requires columns avg, bndr & id
  ## built from the helpful code here https://github.com/darrenkoppel/Metal-toxicity-to-Antarctic-algae/wiki/Dose-response-curves-with-R%2C-plotted-with-ggplot
  df = subset(df, df$id==id)
  model = drm(y~x, data=df, fct=MM.2())
  newdata = expand.grid(conc=exp(seq(log(startval), log(endval), length=150)))
  pm = predict(model, newdata=newdata, interval="prediction")
  newdata$p <- pm[,1]
  newdata$pmin <- pm[,2]
  newdata$pmax <- pm[,3]
  newdata$km = coef(model)[2]
  newdata$error = coef(summary(model))[2,2]
  ## if we label this id column the same as the df id column, ggplot gets very confused and annoyed
  newdata$n_id = id
  return(newdata)
}


conc_mem2particle = function(concs, head_group_area=6.7E-19, hydro_diam=200, surface_area=NULL, leaflets = 2){
  ## give concs in mM!
  ## give diam in nm
  ## gives answer in nM
  if(is.null(surface_area)){
    surface_area = 4*pi*((hydro_diam/1000000000)/2)^2
  }
  out = (((head_group_area)*(concs/1000))/(leaflets*surface_area))*1000000000
  return(out)
}




densPlot = function(df, mf, colours, xlab, ylab, startval=0.001, logax="x"){
  
  plot = ggplot(NULL)+
    geom_point(data=df, aes(colour=id, x = x, y = y))+
    geom_errorbar(data=df, aes(x=x, ymin = y-stdev, ymax=y+stdev, colour=id),  alpha = 0.75, width=0.05)+
    geom_ribbon(data=mf, aes(x=conc, y=p, ymin=pmin, ymax=pmax, fill=n_id), alpha=0.15, colour=NA) +
    geom_line(data=mf, aes(x=conc, y=p, colour=n_id)) +
    theme_bw(base_size =txt_size)+
    scale_colour_manual(values = colours)+
    scale_fill_manual(values=colours)+
    theme(legend.title=element_blank())+
    guides(color = guide_legend(override.aes = list(size=10,pch=20)),
           fill = FALSE)+
    xlab(xlab)+ylab(ylab)+My_Theme
  
  if (tolower(logax)=="x"){
    plot = plot+scale_x_continuous(trans = log_trans(10), labels = comma_format(big.mark = "",
                                                                    decimal.mark = "."))+
      annotation_logticks(sides="b")
  } 
  if (tolower(logax)=="y"){
    plot = plot+scale_y_continuous(trans = log_trans(10), labels = comma_format(big.mark = "",
                                                                                decimal.mark = "."))+
      annotation_logticks(sides="l")
  } 
  #plot <- cowplot::ggdraw(plot) + 
   # theme(plot.background = element_rect(fill="white", color = NA))
  return(plot)
}




















