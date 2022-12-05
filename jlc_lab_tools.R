## jlc lab tools


library(ggplot2)
library(grid)
library(gridExtra)
library(cowplot)



txt_size = 22
lab_size = txt_size*1.222

## theme for plots
My_Theme = theme(
  axis.title.x = element_text(size = txt_size,family="Sans", colour="black"),
  axis.text.x = element_text(size = txt_size,family="Sans", colour="black"),
  axis.title.y = element_text(size = txt_size,family="Sans", colour="black"),
  axis.text.y = element_text(size = txt_size,family="Sans", color="black"),
  strip.text.x = element_text(size = txt_size,family="Sans"),
  strip.text = element_text(size = txt_size,family="Sans"),
  strip.text.y = element_text(size = txt_size,family="Sans"),
  legend.title = element_text(size = txt_size,family="Sans"),
  legend.text = element_text(size = txt_size,family="Sans"),
  strip.background =element_rect(fill="white",color="black"),
  plot.title = element_text(size=txt_size,family="Sans"),
  plot.margin=unit(c(0.6,0.6,0.6,0.6), "cm"),
  panel.background = element_rect(fill = "white", colour = "white"))



norm = function(col){
  ## col    (scalar)    any list of ints to be internally normalised
  col = as.numeric(as.character(col))
  col = (col - min(col))/(max(col)-min(col))
  return(col)
}


avg = function(col_list, name){
  ## col_list   (list)    list of same length columns to have mean average and stdev calculated
  ## name       (str)     identifier to be added to output dataframe
  df = data.frame(rowMeans(data.frame(col_list)),
                  apply(data.frame(col_list),1,sd),
                  name)
  names(df) = c("y", "stdev", "id")
  return(df)
}


xy_data_importr = function(infile, dir, ext, nrml="none", precut=0, columns = c(1:2)){
  ## for importing any single xy dataset
  
  ##  infile    (str)         name of file to import
  ##  dir       (str)         path of directory file is in
  ##  ext       (str)         extension of file eg ".csv"
  ##  nrml      (str)         pass "yes" for full 0-1 normalisation, pass "zero" for zero normalisation
  ##  precut    (int)         remove data before this point
  ##  columns   (lsit, strs)  list of columns to import, defaults to xy, for DLS give c(2,3,5,6,8,9))
  
  ## import data
  s1 <- read.csv(paste0(dir,infile,ext), header=TRUE)
  ## cut to desired columns
  s1 = s1[,columns]
  ## apply names to columns
  nms = rep(c("x", "y"), length(columns)/2)
  names(s1) = nms
  ## remove data before certain point if desired
  s1 = subset(s1, s1$x > precut)
  ## zero x axis
  s1$x = s1$x - s1$x[1]
  
  if (nrml=="zero" | nrml==0){
    s1$y = s1$y - s1$y[1]
  }
  if (nrml=="yes" | nrml=="y"){
    s1$y = norm(s1)
  }
  s1$id = infile
  s1 = na.omit(s1)
  return(s1)
}

data_framr = function(file_list, dir, ext, nrml="n", precut=0){
  ## wrapper for xy_importer, feed in list of files
  df = do.call(rbind, lapply(file_list, xy_data_importr, dir=dir, ext=ext, nrml=nrml, precut=precut))
  return(df)
}


cuvette_importer = function(infile, dir, ext="", nrml="no"){
  ## for importing runs where using multiple cuvettes eg flux data
  
  ##  infile    (str)         name of file to import
  ##  dir       (str)         path of directory file is in
  ##  ext       (str)         extension of file eg ".csv"
  ##  nrml      (str)         pass "yes" for full 0-1 normalisation
  
  
  ## read.csv seems to fail if there are excess blank columns in the csv. maybe worth changing to read.something else?
  s1 = read.csv(paste0(dir, infile, ext), header=TRUE )
  ## remove label lines
  s1 = tail(s1,-1)
  ## drop NA column sometimes added to end
  s1=s1[,!apply(is.na(s1), 2, any)]
  ## make numeric
  s1[] = apply(s1[], 2, function(x) as.numeric(as.character(x)))
  ## remove crud
  s1 = na.omit(s1)
  ## make vector of cuvette labels, max 4
  samples = na.omit(names(s1))
  samples = samples[seq(1,length(samples),2)]
  ## convert to longform
  outdf = data.frame(x=NA,y=NA,sample=NA)
  for (i in seq(1,(length(samples)*2-1), 2)){
    
    ## cuts data to a single xy column, and normalises y column if requested
    tmp = data.frame(x=s1[,i], y=s1[,i+1], sample=samples[[(i+1)/2]])
    if (nrml=="y" | nrml=="yes"){
      tmp[,2]=norm(tmp[,2])*100
    }
    outdf = rbind(outdf, tmp)
  }
  print("ignore the errors and check to see if data is in your enviornment now")
  outdf = na.omit(outdf)
  outdf$id = infile
  return(outdf)
}


cuvette_framer = function(file_list, dir, ext="", nrml="no"){
  ## wrapper for cuvette importer to import multiple files
  df = do.call(rbind, lapply(file_list, cuvette_importer, dir=dir, ext=ext, nrml=nrml))
  return(df)
}


plotr = function(df, xlab="x", ylab="y", leg_loc="right", error="n",
                 cols=NULL, fac=NULL, marker="line", ticks=10, log="no"){
  # makes scatter plots
  
  ##    df      (dataframe)   data to plot
  ##    xlab    (str)         x axis label
  ##    ylab    (str)         y axis label
  ##    leg_loc (str)         location of legend relative to plot eg "left" or "top"
  ##    error   (str)         give "y" to add errorbars to plot
  ##    cols    (list, str)   colours to plot in, same order as labels
  ##    fac     (str)         column title to facet plots with eg. "id"
  ##    marker  (str)         line or point 
  ##    ticks   (int)         number of axis plots
  ##    log     (str)         feed "y" or "x" to log axis

  
  if (error == "yes" | error == "y"){
    plot = ggplot(df, aes(x=x, y=y, colour=id,fill=id,
                          ymax=y+stdev, ymin=y-stdev))+
      geom_ribbon(alpha = 0.15, colour=NA)
    
  }else if (error=="half"|error==0.5){
    plot = ggplot(df, aes(x=x, y=y, colour=id,fill=id,
                          ymax=y+(0.5*stdev), ymin=y-(0.5*stdev)))+
      geom_ribbon(alpha = 0.15, colour=NA)
    
  } else{
    plot = ggplot(df, aes(x=x, y=y, colour=id,fill=id))
      
  }
  
  plot = plot + theme_bw(base_size = txt_size)+
    theme(text = element_text(size=txt_size))+
    xlab(xlab)+ylab(ylab)+
    My_Theme+
    theme(legend.position=leg_loc,
          legend.title=element_blank())+
    scale_x_continuous(breaks = scales::pretty_breaks(n = ticks))+ 
    guides(color = guide_legend(override.aes = list(size=10,pch=20)),
           fill = FALSE)+
    
    geom_blank()
  
    if (marker=="line"){
    plot = plot + 
      geom_path()
    }else{
      plot = plot+
       geom_point()
    }
    if(is.null(cols)){
      plot = plot
    }else{
      plot = plot +    
        scale_colour_manual(values = (cols))+
        scale_fill_manual(values = (cols))
    }
    if(is.null(fac)){
      plot=plot
    }else{
      plot = plot+facet_grid(get(fac)~.)
    }
    if (log=="x"){
      library(scales)
      plot = plot+
        scale_x_continuous(trans = log_trans(10), labels = comma_format(big.mark = "",
                                                                        decimal.mark = "."))+#,limits = log_lim, breaks=c(1,10,100,1000,1000,10000,10000)) +
        annotation_logticks(sides="b")
    } else if (log=="y"){
      library(scales)
      plot = plot+
        scale_y_continuous(trans = log_trans(10), labels = comma_format(big.mark = "",
                                                                        decimal.mark = "."))+#,limits = log_lim, breaks=c(1,10,100,1000,1000,10000,10000)) +
        annotation_logticks(sides="b")
    }  
    return(plot)
}


barPlotr = function(bardf, labels, cols, ylab, error){
  # makes bar plots
  
  ##    bardf   (dataframe)   data to plot
  ##    cols    (list, str)   colours to plot in, same order as labels
  ##    ylab    (str)         y axis label
  ##    error   (str)         give "y" to add errorbars to plot
  if (error == "yes" | error == "y"){
    plot = ggplot(bardf, aes(x=id, y=y, colour=id,fill=id,
                          ymax=y+stdev, ymin=y))+
      geom_errorbar(position="dodge", width=0.2, colour="black")
      
    
  }else if (error=="half"|error==0.5){
    plot = ggplot(bardf, aes(x=id, y=y, colour=id,fill=id,
                          ymax=y+(0.5*stdev), ymin=y))+
      geom_errorbar(position="dodge", width=0.2, colour="black")
      
  } else{
    plot = ggplot(bardf, aes(x=id, y=y, colour=id,fill=id))
  }
    plot = plot + geom_bar(stat="identity",
                           position=position_dodge(), alpha = 0.1, size=1.5) +
      scale_fill_manual(values=cols,
                        labels  = labels)+
      scale_colour_manual(values=cols)+
      guides(fill = FALSE)        +
      guides(color = FALSE)        +
      theme(legend.position="none")+
      theme_bw(base_size=txt_size)+
      xlab("")+ylab(ylab)+
      My_Theme

  return(plot)
}







effluxDataGrab = function(df, labels, injects, wavelengths=1, nrml="y"){
  # prepares imported flux data for plotting
  
  ##  df          (dataframe)     output from cuvette_framer tool
  ##  labels      (list, str)     labels for each sample, need to factor onto id column before running this
  ##  injects     (list, int)     injection points for each sample, in same order as labels
  ##  wavelengths (int)           set to 2 if doing dual wavelength work eg. HTPS flux analyis
  ##  nrml      (str)         pass "yes" for full 0-1 normalisation
  

  outdf = data.frame()
  for (i in seq(1, length(labels))){
    ## subset to single run on the fluorimeter
    tmp = subset(df, df$id==labels[[i]])
    col_list = list()
    x_list = list()
    
    if (wavelengths == 2){
      for (j in seq(1, length(unique(tmp$sample)),2)){
        ## subset to each cuvette and each wavelength, then subtract one 
        ## from another for corrections
        t2 = subset(tmp, tmp$sample==unique(tmp$sample)[[j]])  
        t3 = subset(tmp, tmp$sample==unique(tmp$sample)[[j+1]])
        ## remove pre-inject data prior to normalisation
        t2 = subset(t2, t2$x > injects[[i]])
        t3 = subset(t3, t3$x > injects[[i]])
        ## subtract one wavelength from t'other
        t2y = t3$y - t2$y
        ## normalise
        if (nrml=="yes" | nrml=="y"){
          t2y = norm(t2y)*100
        }
        ## add to list for averaging
        col_list[[j]]=t2y
        x_list[[j]]=t2$x

               
      }
    } else {    
    for (j in seq(1, length(unique(tmp$sample)))){
      ## subset to each cuvette
      t2 = subset(tmp, tmp$sample==unique(tmp$sample)[[j]])
      ## remove pre-inject data prior to normalisation
      t2 = subset(t2, t2$x > injects[[i]])
      t2y = t2$y
      ## normalise
      if (nrml=="yes" | nrml=="y"){
        t2y = norm(t2y)*100
      }
      col_list[[j]]=t2y
      x_list[[j]]=t2$x    
      }
    }
    ## remove empty entries due to counting in 2's in  loop
    col_list <- col_list[!sapply(col_list,is.null)]
    x_list <- x_list[!sapply(x_list,is.null)]
    ## calculate average and sd
    a1 = try(avg(col_list, labels[[i]]))
    ## this part of the function will step in if the different runs are different lengths (e.g have been re-cut)
    ## because they are already normalised, we should only need to cut each one down to smallest length
    if (class(a1) == "try-error"){
      cut = min(lengths(col_list))
      col_list = lapply(col_list, function(x) x=x[1:cut])
      x_list = lapply(x_list, function(x) x=x[1:cut])
      a1 = (avg(col_list, labels[[i]]))
      print("ignore the error everything is fine")
      print("if it's not fine, have you remembered to factor the labels column from your input df?")
    }
    ## add time data back in
    time = rowMeans(data.frame(x_list))
    a1$x = time
    ## zero normalise & ints
    a1$x = a1$x-a1$x[1]
    outdf = rbind(a1, outdf)
  }
  return(outdf)
}

maxFlux = function(df, labels, tmax){
  ## calculates background corrected maximum values fot efflux data
  
  ##  df          (dataframe)     output from effluxDataGrab tool
  ##  labels      (list, str)     labels for each sample, need to factor onto id column before running this
  # labels should be in order bkg1, sample1, bkg2, sample2, etc etc
  ##  tmax        (int)           timepoint to calculate max values at
  indf = subset(df, df$x < tmax)
  bardf = data.frame(y=indf$y, stdev=indf$stdev, id=indf$id)
  bkg1 = subset(bardf, bardf$id == labels[1])
  bkg1 = subset(bkg1, bkg1$y == max(bkg1$y))
  dat1 = subset(bardf, bardf$id == labels[2])
  dat1 = subset(dat1, dat1$y == max(dat1$y))
  bkg2 = subset(bardf, bardf$id == labels[3])
  bkg2 = subset(bkg2, bkg2$y == max(bkg2$y))
  dat2 = subset(bardf, bardf$id == labels[4])
  dat2 = subset(dat2, dat2$y == max(dat2$y))
  ## stdev calculation
  outdf = rbind(data.frame(id=labels[2], y=(dat1$y-bkg1$y),
                           stdev = sqrt(dat1$stdev+bkg1$stdev)),
                data.frame(id=labels[4], y=(dat2$y-bkg2$y),
                           stdev = sqrt(dat2$stdev+bkg2$stdev)))

  return(outdf)
}

initRate = function(df, rates=c(1,1.5,2.5)){
  ## for calculating initial rates
  ##  df          (dataframe)     output from effluxDataGrab tool
  ##  rates       (list, ints)    3x time values that data will be reduced to and linear fit to. Then averaged to give initial rate
  init_df=data.frame()
  for (i in unique(df$id)){
    t1 = subset(df, df$id==i)
    r1 = lm(y~x, subset(t1, t1$x < rates[1]))
    r2 = lm(y~x, subset(t1, t1$x < rates[2]))
    r3 = lm(y~x, subset(t1, t1$x < rates[3]))  
    tmp = data.frame(y=mean(c(r1$coefficients[2],r2$coefficients[2],r3$coefficients[2])),
                     stdev=sd(c(r1$coefficients[2],r2$coefficients[2],r3$coefficients[2])),
                     id=i)
    init_df = rbind(init_df, tmp)
  }
  return(init_df)
}

## not my tool but very useful
find_peaks <- function (x, m = 3){
  ## x  (scalar)     vector you wish to find peaks of
  ## m  (int)       sensitivity, larger m = fewer peaks
  # for valley detection, pass -x.
  shape <- diff(sign(diff(x, na.pad = FALSE)))
  pks <- sapply(which(shape < 0), FUN = function(i){
    z <- i - m + 1
    z <- ifelse(z > 0, z, 1)
    w <- i + m + 1
    w <- ifelse(w < length(x), w, length(x))
    if(all(x[c(z : i, (i + 2) : w)] <= x[i + 1])) return(i + 1) else return(numeric(0))
  })
  pks <- unlist(pks)
  pks
}

                      pngOut <- function(plot, dir, name, width=210, height=180){
  ##  plot    (obj)     graphic to save to png
  ##  dir     (str)     abs path to save png to
  ##  name    (str)     filename to save plot as
  ##  width   (int)     width of final output (mm)
  ##  height   (int)    height of final output (mm)
  
  ## apply white background so no werid box appears 
  out <- cowplot::ggdraw(plot) + 
    theme(plot.background = element_rect(fill="white", color = NA))
  loc = paste0(dir, name, ".png")
  ## save
  ggsave(filename=loc, 
         plot = out, 
         device = png, 
         width = width, 
         height = height, 
         units = "mm")
}
  

deriviative1 <- function(x,y){
  ## Performs 1st derivative calculation on XY data, returns data frame with original x axis and 1st diff Y
  ## x    (col, int)    x axis column to differentiate
  ## y    (col, int)    y axis colun to differentiate
  tdf = data.frame(x,y)
  model = smooth.spline(x=tdf[,1], y=tdf[,2])
  ddf = data.frame(predict(model, x=x, deriv=1))
  return(ddf[,2])
}
