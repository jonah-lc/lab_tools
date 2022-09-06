#### Functions for plotting densitometry of gels #####se
library(ggplot2)
library(drc)
library(scales)

source("/media/jonah/Lickitung/SYNC/4-SCRIPTS/R/jlc_lab_tools.R")

txt_size = 22



idensImport = function(path, id) {
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
  ## imported_df    (df)    output from "import" module
  
  ## background correction
  imported_df$corrected = imported_df$ints - imported_df$bkg
  ## normalise
  imported_df$norm = norm(imported_df$corrected)
  ## turn to percentage
  imported_df$pb = (1-imported_df$norm)*100
  return(imported_df)
}



preAvg = function(folder, file, ext=".csv", binder_conc){
  inp = import(paste0(folder, file, ext))
  inp = perc_bound(inp)
  inp$binder = binder_conc
  ## now remove last line post normalisation
  inp = inp[1:(nrow(inp)-1), ]
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



densCalc= function(folder, file_list, ext=".csv", binder_conc){
  # binder_conc     (scalar)      list of int for concentration of titrator, must be same half length of data files
  
  ## prepare data
  df = do.call(rbind, lapply(file_list, preAvg, folder=folder, ext=ext, binder_conc))
  # how many repeats
  dat_numb = length(unique(df$id))
  ## avg and sd
  avg = aggregate(pb ~ binder, df, mean)
  stdev = aggregate(pb ~ binder, df, sd)
  ## put back into dataframe
  df$avg = rep(rep(avg[,2], dat_numb))
  df$stdev = rep(rep(stdev[,2], dat_numb))
  ## make common id for avgd data
  df$c_id = comSub(df$id)
  return(df)
}
  
fitter = function(in_id, df, length_out=1000){
  ## will fit model to data from dens_cals.
  ## df requires columns avg, bndr & id
  df = subset(df, df$c_id==in_id)
  tmp = data.frame(avg = df$avg, bndr=df$binder, id=in_id)
  tmp = unique(tmp)
  model.drm <- drm (avg ~ bndr, data = tmp, fct = MM.2())
  mml <- data.frame(S = seq(0, max(tmp$bndr), length.out = length_out))
  print(summary(mml))
  mml$v <- predict(model.drm, newdata = mml)
  mml$c_id = unique(tmp$id)
  return(mml)
}

fitCalc = function(df, id_list){
  # df should be output of densCalc
  # id_list can be unique(df$c_id)
  out = do.call(rbind, lapply(id_list, fitter, df))
  return(out)
}

kdr = function(df){
  ## df = output from densCals
  kdf = data.frame()
  for (h in unique(df$c_id)){
    kdlist = list()
    tf = subset(df, df$c_id==h)
    for (i in seq(1, length(unique(tf$id)))){
      tid = unique(tf$id)[[i]]
      tmp = subset(df, df$id==tid)
      mf1 = fitter(in_id = tid, df = data.frame(avg=tmp$pb, binder=tmp$binder, c_id=tmp$id))
      kd = mean(mf1$S[round(mf1$v,0)==50])
      kdlist = append(kdlist, kd)
    }
    
    kd = mean(sapply(kdlist, mean))
    kdstdev = sd(sapply(kdlist, mean))
    out = data.frame(kd=kd, stdev=kdstdev, id=h)
    kdf = rbind(kdf, out)
    }
  return(kdf)
}

conc_mem2particle = function(concs, head_group_area=6.7E-19,hydro_diam=200, leaflets=2){
  ## give concs in mM!
  ## give diam in nm
  ## gives answer in nM
  surface_area = 4*pi*((hydro_diam/1000000000)/2)^2
  out = (((head_group_area)*(concs/1000))/(leaflets*surface_area))*1000000000
  return(out)
}

densPlot = function(df, mf, labs, cols, xlab, ylab, log="no", colvar="df$c_id"){
  ## colvar = df column name to colour by 
  # correct labels
  df$c_id = factor(df$c_id, levels = unique(df$c_id), labels = labs)
  mf$c_id = factor(mf$c_id, levels = unique(mf$c_id), labels = labs)
  ## set colvar
  mvar = paste0("mf",substr(colvar, 3, nchar(colvar)))
  # prepare plot
  out = ggplot()+
    geom_point(data=df, aes(x=df$binder, y=df$avg, color = df$var2),   size = 1)+
    geom_errorbar(data=df, aes(x=df$binder, ymin=df$avg-df$stdev, ymax=df$avg+df$stdev,
                                color = df$var2),  alpha = 0.75, width=0.1)+
    geom_line(data = mf, aes(x=mf$S, y=mf$v, colour = mf$var2) , size = 1)+
    theme_bw(base_size =txt_size)+
    xlab(xlab)+ylab(ylab)+
    theme(legend.position="right",
         strip.background = element_blank())+
    scale_colour_manual(values = cols, name="")
    

    if (log=="yes"|log=="y"){
      out = out + scale_x_continuous(trans='log10', 
                                     breaks=trans_breaks('log10', function(x) 10^x),
                                                        labels=trans_format('log10',
                                                                            math_format(10^.x)))
    }  
    
    return(out)
}






