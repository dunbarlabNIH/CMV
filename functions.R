gfp_plot=function(...,cell_type=1,monkeys){
  x=list(...)
  nmonkey=length(x)
  
  df=data.frame(time=NULL,GFP=NULL,Monkey=NULL,color=NULL)
  for(i in 1:nmonkey){
    temp=data.frame(time=rownames(x[[i]]),GFP=x[[i]][,cell_type],Monkey=monkeys[i],color=monkeys[i])
    if(monkeys[i]=="jm82"){
      temp$color="blue"
    }else if(monkeys[i]=="jd76"){
      temp$color="red"
    }else{
      temp$color="black"
    }
    df=rbind(df,temp)
  }
  
  df=df[!is.na(df$GFP),]
  print(ggplot(df,aes(x=as.numeric(as.character(time)),y=GFP,color=Monkey,group=Monkey))+
          geom_line(color=df$color)+
          geom_point(color=df$color)+scale_x_continuous("Months post-transplant",limits=c(0,20))+scale_y_continuous("GFP Percent")+theme_grey(base_size=22))
}

name <- function(v1) {
  deparse(substitute(v1))
}


barcode_ggheatmapN_bar <- function(your_barcoding_data,
                              your_facs_data,
                              names = colnames(your_barcoding_data),
                              n_clones = 10,
                              your_title = "",
                              grid = TRUE,
                              label_size = 1,
                              dendro = FALSE,
                              cellnote_size = 4,
                              printtable = FALSE,
                              table_option = "percents",
                              log_transform = TRUE,
                              log_choice = exp(1),
                              distance_method = "Euclidean",
                              minkowski_power = 1,
                              cellnote_option = "stars",
                              hclust_linkage = "complete",
                              row_order = "hierarchical",
                              clusters = 0,input=NA) {
  your_barcoding_data_list <- list(raw_reads = your_barcoding_data, prop_table = as.data.frame(prop.table(as.matrix(your_barcoding_data),2)))
  your_barcoding_data_list$prop_table=t(t(your_barcoding_data_list$prop_table)*your_facs_data)
  your_barcoding_data_list$prop_table[is.na(your_barcoding_data)] <- 0
  max=max(your_barcoding_data_list$prop_table)
  if(log_transform){
    max=log(max,base=log_choice)
  }
  if (any(colSums(your_barcoding_data_list$prop_table) == 0)){
    stop("One of your columns contained no data")
  }
  your_barcoding_data_list$prop_table[your_barcoding_data_list$prop_table == 0] <- NA
  #creates data frame that shows rank of original your_barcoding_data
  your_barcoding_data_list <- c(your_barcoding_data_list, list(ranks = apply(-your_barcoding_data_list$prop_table, 2, rank, ties.method = "min", na.last = "keep")))
  #subsets those barcodes that have at least one top N clone
  top_clones_choices <- apply(your_barcoding_data_list$ranks, 1, function(x){any(x<=n_clones, na.rm = TRUE)})
  your_barcoding_data_list <- lapply(your_barcoding_data_list, function(x) {x[top_clones_choices,]})
  #creates data frame with '*' for those cells w/ top clones
  your_barcoding_data_list <- c(your_barcoding_data_list, list(cellnote_matrix = your_barcoding_data_list$ranks))
  your_barcoding_data_list$cellnote_matrix[your_barcoding_data_list$cellnote_matrix > n_clones] <- NA
  your_barcoding_data_list$cellnote_matrix[your_barcoding_data_list$cellnote_matrix <= n_clones] <- "*"
  your_barcoding_data_list$prop_table[is.na(your_barcoding_data_list$prop_table)] <- 0
  #takes log of data
  your_barcoding_data_list <- c(your_barcoding_data_list, list(logged_data = custom_log(your_barcoding_data_list$prop_table, log_choice, FALSE)))
  if(row_order == "hierarchical") {
    if(distance_method == "Minkowski"){
      hclustering <- hclust(proxy::dist((if (log_transform) your_barcoding_data_list$logged_data else your_barcoding_data_list$prop_table), method = distance_method, p = minkowski_power), method = hclust_linkage)
    } else {
      hclustering <- hclust(proxy::dist((if (log_transform) your_barcoding_data_list$logged_data else your_barcoding_data_list$prop_table), method = distance_method), method = hclust_linkage)
    }
    barcode_order <- rev(hclustering$order)
  } else if(row_order == "emergence") {
    barcode_order <- do.call(order, as.data.frame(your_barcoding_data_list$prop_table))
  } else {
    stop("row_order must be one of \" hierarchical\" or \"emergence\"")
  }
  
  if(log_transform){
    plotting_data <- reshape2::melt(as.matrix(your_barcoding_data_list$logged_data))
    actual_scale <- c(log(100/4000000, log_choice) - 1,mean(c(log(100/4000000, log_choice) - 1,max)) ,max)
    your_scale <- scales::rescale(actual_scale, to = c(0,1))
    your_breaks <- c(actual_scale)
    your_labels <- c("Defined 0","", "")
    your_limits <- c(log(100/4000000, log_choice)-1,max)
    your_colors <- c("#white", "#gray", "black")
  } else {
    plotting_data <- reshape2::melt(as.matrix(your_barcoding_data_list$prop_table))
    your_scale <- c(0, max(your_barcoding_data_list$prop)/2,max(your_barcoding_data_list$prop))
    your_breaks <- c(0, .5,1)
    your_labels <- c("0%","", "max")
    your_limits <- c(0,1)
    your_colors <- c("#4575B4", "#fefeb9", "red4")
  }
  
  
  
  colnames(plotting_data) <- c("BARCODE", "SAMPLE", "Size")
  plotting_data$BARCODE <- factor(plotting_data$BARCODE, levels = rev(rownames(your_barcoding_data_list$prop_table)[barcode_order]))
  input=as.vector(input[unique(plotting_data$BARCODE),])
  print(input)

  plotting_data=rbind(data.frame(BARCODE=length(unique(plotting_data$BARCODE)),SAMPLE="NA",Size=NA),plotting_data)
  plotting_data$color=c(input,rep("",length(input)*(length(names)-2)+1))
  if (cellnote_option %in% c("logs", "ranks", "stars", "percents", "reads")){
    if(cellnote_option == "logs"){
      plotting_data$CELLNOTE <- round(reshape2::melt(as.matrix(your_barcoding_data_list$logged_data))$value, 2)
      cellnote_vjust = 0.5
    }
    if(cellnote_option == "ranks"){
      plotting_data$CELLNOTE <- reshape2::melt(as.matrix(your_barcoding_data_list$ranks))$value
      cellnote_vjust = 0.5
    }
    if(cellnote_option == "stars"){
      plotting_data$CELLNOTE <- c(rep("*",length(input)),rep(NA,nrow(plotting_data)-length(input)))
      cellnote_vjust = 0.75
    }
    if(cellnote_option == "percents"){
      plotting_data$CELLNOTE <- paste0(round(reshape2::melt(as.matrix(your_barcoding_data_list$prop_table))$value*100,2), "%")
      cellnote_vjust = 0.5
    }
    if(cellnote_option == "reads"){
      plotting_data$CELLNOTE <- reshape2::melt(as.matrix(your_barcoding_data_list$raw_reads))$value
      cellnote_vjust = 0.5
    }
  } else {
    stop("cellnote_option must be one of c(\"logs\", \"ranks\", \"stars\", \"percents\", \"reads\")")
  }
  
  if(grid)gridColor = "black" else gridColor = NA
  
  
  false_column_label <- unique(colnames(your_barcoding_data_list$logged_data))[which(max(nchar(unique(colnames(your_barcoding_data_list$logged_data)))) == nchar(unique(colnames(your_barcoding_data_list$logged_data))))[1]]
  
  if(printtable == TRUE){
    switch(table_option,
           reads = return(your_barcoding_data_list$raw_reads[barcode_order,]),
           percents = return(your_barcoding_data_list$prop_table[barcode_order,]),
           logs = return(your_barcoding_data_list$logged_data[barcode_order,]),
           ranks = return(your_barcoding_data_list$ranks[barcode_order,]))
  } else {
    
    g1_heatmap <- ggplot2::ggplot(plotting_data, ggplot2::aes(x = SAMPLE, y = BARCODE,colour=Size))+
      ggplot2::geom_tile(ggplot2::aes(fill = Size),colour=NA)+
      ggplot2::geom_text(ggplot2::aes(label = CELLNOTE), vjust = cellnote_vjust, size = cellnote_size)+
      ggplot2::scale_fill_gradient(low="white",high="purple4")+
      ggplot2::scale_y_discrete(labels = NULL, breaks = NULL, expand = c(0,0))+
      ggplot2::scale_x_discrete(expand = c(0,0), labels = names)+
      ggplot2::ylab(NULL)+
      ggplot2::xlab(NULL)+
      ggplot2::ggtitle(paste0("\n", your_title, "\n"))+
      ggplot2::theme(
        plot.title = ggplot2::element_text(size = 20),
        axis.text.x = ggplot2::element_text(angle=90, hjust = 1, vjust = 0.5, size = label_size),
        legend.text = ggplot2::element_text(size =  15, face = 'bold'),
        legend.title = ggplot2::element_text(size =  15),
        axis.ticks = ggplot2::element_blank(),
        legend.key.height=ggplot2::unit(5,"line"))
    
    
    if(row_order == 'emergence'){
      print(g1_heatmap)
    } else {
      dendro_data <- ggdendro::dendro_data(hclustering, type = 'rectangle')
      g2_dendrogram <- ggplot2::ggplot(ggdendro::segment(dendro_data))+
        ggplot2::geom_segment(ggplot2::aes(x = x, y = y, xend = xend, yend = yend))+
        ggplot2::scale_x_discrete(expand = c(.5/nrow(your_barcoding_data_list$logged_data),0.01))+
        ggplot2::scale_y_reverse(expand = c(0.01,0), labels = false_column_label, breaks = if(row_order == 'emergence') NULL else mean(hclustering$height))+
        ggplot2::coord_flip()+
        ggplot2::ylab(NULL)+
        ggplot2::xlab(NULL)+
        ggplot2::ggtitle("\n \n")+
        ggplot2::theme(
          plot.title = ggplot2::element_text(size = 20),
          axis.text.x = ggplot2::element_text(colour = 'white', angle=90, hjust = 1, vjust = 0.5, size = label_size),
          panel.background = ggplot2::element_rect(fill = "white",colour = "white"),
          axis.ticks = ggplot2::element_blank()
        )
      
      if(dendro){
        
        if(clusters > 0){
          # custom_colors <- c(RColorBrewer::brewer.pal(8, 'Set1'), "cyan", "black", "grey")
          custom_colors <- c("#89C5DA", "#DA5724", "#74D944", "#CE50CA", "#3F4921", "#C0717C", "blue")
          clusters_df <-data.frame(CLUSTERS = cutree(hclustering, clusters), ASSIGNMENT=factor(hclustering$labels,levels=hclustering$labels[(hclustering$order)]))
          g3_clusters <-ggplot2::ggplot(clusters_df, ggplot2::aes(x =1,y=ASSIGNMENT,fill=factor(CLUSTERS)))+
            ggplot2::geom_tile()+
            ggplot2::ggtitle("\n \n")+
            ggplot2::scale_fill_manual(values = custom_colors[1:clusters])+
            ggplot2::scale_x_continuous(expand=c(0,0), labels = false_column_label, breaks = 1)+
            ggplot2::theme(
              plot.title = ggplot2::element_text(size = 20),
              axis.title=ggplot2::element_blank(),
              axis.ticks=ggplot2::element_blank(),
              axis.text.y=ggplot2::element_blank(),
              axis.text.x = ggplot2::element_text(colour = 'white', angle=90, hjust = 1, vjust = 0.5, size = label_size),
              legend.position="none")
          gridExtra::grid.arrange(g2_dendrogram,
                                  g3_clusters,
                                  g1_heatmap,
                                  ncol = 3,
                                  widths = c(1,.2,4))
          
        } else {
          gridExtra::grid.arrange(g2_dendrogram,
                                  g1_heatmap,
                                  ncol = 2,
                                  widths = c(1,4))
        }
        
        
      } else {
        g1_heatmap
      }
    }
    
    
    
    
  }
  
}


custom_log <- function(x, log_choice, vary_log_min){
  x <- log(x, log_choice)
  #sets the -inf values to be the minimum of the possible value in the data frame -1
  if(vary_log_min) {
    x[x == -Inf] <- (min(x[x > -Inf]) - 1)
  } else {
    x[x == -Inf] <- (log(100/4000000, log_choice)-1)
  }
  return(x)
}

cbc=function(cbc,day,legend=FALSE){
  ntime=nrow(cbc)
  df=NULL
  
  df$pbmnc=cbc$WBC*(cbc$LYMP+cbc$MONOS)/100
  df$monos=cbc$WBC*cbc$MONOS/100
  df$segs=cbc$WBC*cbc$SEGS/100
  df$lymp=cbc$WBC*cbc$LYMP/100
  df$cd3=cbc$CD3/cbc$P2*df$pbmnc
  df$cd8=cbc$CD8/cbc$P2*df$pbmnc

  df$cd4=cbc$CD4/cbc$P2*df$pbmnc
  df$cd20=cbc$CD20/cbc$P2*df$pbmnc

  df$nkg2a=cbc$NKG2A/cbc$P2*df$pbmnc
  df$cd16=cbc$CD16/cbc$P2*df$pbmnc
  df$cd56=cbc$CD56/cbc$P2*df$pbmnc
  df$rbc=cbc$RBC
  df$plt=cbc$PLT
  df$wbc=cbc$WBC
  df$grans=cbc$SEGS*cbc$WBC/100
  
  df$time=cbc$X.Daysinfus-day
  df=as.data.frame(df)
  
  cbc.df=melt(df[,c("cd4","cd8","cd16","cd56")])
  cbc.df$time=rep(df$time,times=4)
  
  cbc.df$value=cbc.df$value*1000
  print(ggplot(cbc.df,aes(time,value,color=variable))+
    geom_point(size=3)+
    geom_line(size=2)+scale_x_continuous("Days post rhCMV infection",breaks=seq(-30,1000,by=30),limits=c(-30,225))+
    scale_y_continuous("Number of cells per uL PB",trans = log10_trans(),breaks = trans_breaks("log10", function(x) 10^x),limits=c(1,10^6),labels = trans_format("log10", math_format(10^.x)))+ scale_colour_discrete(name="Cell type")+
    scale_color_manual("Cell Type",labels=c("T CD4+","T CD8+","NK CD56-CD16+","NK CD56+CD16-"),values=c("black","red","green","blue"))+theme_grey(base_size=22)+theme(legend.position=ifelse(legend,"right","none")))
  return(cbc.df)
}

cbc.pre=function(cbc,day){
  ntime=nrow(cbc)
  df=NULL
  
  df$pbmnc=cbc$WBC*(cbc$LYMP+cbc$MONOS)/100
  df$monos=cbc$WBC*cbc$MONOS/100
  df$segs=cbc$WBC*cbc$SEGS/100
  df$lymp=cbc$WBC*cbc$LYMP/100
  df$cd3=cbc$CD3/cbc$P2*df$pbmnc
  df$cd20=cbc$CD20/cbc$P2*df$pbmnc
  
  df$cd16=cbc$CD16/cbc$P2*df$pbmnc
  df$cd56=cbc$CD56/cbc$P2*df$pbmnc
  df$rbc=cbc$RBC
  df$plt=cbc$PLT
  df$wbc=cbc$WBC
  df$grans=cbc$SEGS*cbc$WBC/100
  
  df$time=cbc$X.Daysinfus-day
  df=as.data.frame(df)
  df
  cbc.df=melt(df[,c("cd20","cd3","cd16","cd56")])
  cbc.df$time=round(rep(df$time,times=4)/30)
  
  cbc.df$value=cbc.df$value*1000
  
#  ggplot(cbc.df,aes(time,value,color=variable))+
#    geom_point(size=3)+
#    geom_line(size=2)+scale_x_continuous("Months post transplant",limits=c(0,6),breaks=seq(1,6))+
#    scale_y_continuous("Number of cells per uL PB",trans = log10_trans(),breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x)),limits=c(1,10^4))+ scale_colour_discrete(name="Cell type")+
#    scale_color_manual(labels=c("B","T","NK CD56-CD16+","NK CD56+CD16-"),values=c("black","red","green","blue"))+theme_grey(base_size=22)
  return(cbc.df)
}
barplot.barcode=function(data,columns,tbgrans,months,width){
  
  proptable=apply(data[,c(columns)],2,function(x){x/sum(x)})
  to.compare=tbgrans
  
  biased.barcodes=rownames(proptable)[apply(proptable[,columns]>10*to.compare,1,any)]
  biased.barcodes=biased.barcodes[apply(proptable[biased.barcodes,columns],1,max)>.01]
  
  printtable=proptable[biased.barcodes,columns]
  
  maximum=max(colSums(printtable))
  cols=colSums(printtable)
  barcodes=biased.barcodes
  colnames(printtable)=months
  library(RColorBrewer)
  n <- nrow(printtable)
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  color = color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
  col_vector = c(sample(color,n,replace=T),"white")
  #pie(rep(1,n), col=sample(col_vector, n))
  
  
  colors=as.data.frame(col_vector);colnames(colors)="colors"
  
  
  rownames(colors)=c(barcodes,"other")
  printtable=rbind(printtable,NA)
  
  rownames(printtable)=rownames(colors)
  
  
  printtable["other",]=as.vector(maximum-cols)
  
  
  temp=melt(printtable)
  temp$barcode=rep(unique(rownames(printtable)),ncol(printtable))
  temp$color=rep(colors$colors,ncol(printtable))
  
  
  bar=nrow(printtable)
  barcode=rownames(printtable)
  months=as.numeric(months)
  temp$variable=as.numeric(rep(months,each=nrow(printtable)))
  
  plusminus=width/2
  
  poly=data.frame(x=c(months[1]+plusminus,months[2]-plusminus,months[2]-plusminus,months[1]+plusminus),
                  y=c(maximum,maximum,maximum-temp$value[temp$variable==months[2]&temp$barcode=="other"],maximum-temp$value[temp$variable==months[1]&temp$barcode=="other"]),
                  group="other")
  
  temp=temp[,c("value","barcode","color","variable")]
  
  colnames(temp)=c("Percent","Barcode","Color","Month")
  
  temp2=rbind(temp,temp)
  temp2[1:nrow(temp),"Month"]=temp2$Month[1:nrow(temp)]+width/2
  temp2[(nrow(temp)+1):(2*nrow(temp)),"Month"]=temp2$Month[(nrow(temp)+1):(2*nrow(temp))]-width/2
  
  temp=temp[order(temp$Month),]
  temp2=temp2[order(temp2$Month),]
  
  line=data.frame(x=c(rep(months-width/2,each=2),rep(months+width/2,each=2)),y=rep(c(0,maximum),times=length(months)*4),group=rep(1:(length(months)*2),each=2))
  ggplot()+
    geom_area(data=temp2,mapping=aes(Month,Percent,fill=(Color)),
              stat="identity",position=position_stack(reverse=T))+
    scale_fill_manual(values=as.vector(temp2$Color))+
    #geom_bar(data=temp,mapping=aes(Month,Percent,fill=(Color)),
    #stat="identity",fill=as.vector(temp$Color),width =width)+
    #geom_polygon(data=poly,aes(x=x,y=y,fill=group),fill="white")+
    geom_line(data=line,aes(x=x,y=y,group=group))+
    scale_y_continuous("Percent Contribution",limits=c(0,max(colSums(printtable))),breaks=seq(0,1,by=.1),labels=paste(seq(0,100,by=10),"%"))+
    scale_x_continuous()+
    theme(legend.position="none")
}

tbgrans.max=function(...,months=1:9){
  x=list(...)
  
  for(i in months){
    for(j in 1:length(x)){
      if(!(i %in% colnames(x[[j]]))){
        
        temp.colnames=colnames(x[[j]])
        
        x[[j]]=cbind(x[[j]],NA)
        colnames(x[[j]])=c(temp.colnames,i)
      }
    }
  }
  
  x[[1]]=x[[1]][,order(colnames(x[[1]]))];x[[2]]=x[[2]][,order(colnames(x[[2]]))];x[[3]]=x[[3]][,order(colnames(x[[3]]))]
  max=matrix(ncol=ncol(x[[1]]),nrow=nrow(x[[1]]))
  
  for(i in 1:ncol(max)){
    for(j in 1:nrow(max)){
      max[j,i]=max(x[[1]][j,i],x[[2]][j,i],x[[3]][j,i],na.rm=TRUE)
    }
  }
  
  return(max)
}

cbc.unique=function(cbc){
  cbc=cbc[,c("X.Daysinfus","Monkey")]
  cbc=rev(paste(cbc[,1],cbc[,2]))
  
  return(rev(!duplicated(cbc)))
}

biased_barcode_barplot <- function(your_data, n_timepoints, n_celltypes,timepoints, need_table = FALSE,frac_cutoff=.01,fold_cutoff=10){
  samples_for_plotting <- colnames(your_data)[(ncol(your_data)-n_timepoints+1):ncol(your_data)]
  prop_table <- as.data.frame(prop.table(as.matrix(your_data), margin = 2))
  biased_barcodes_in_each_sample_only <- lapply(1:n_timepoints, function(i) {
    temp <- prop_table[,seq(i, ncol(prop_table), by = n_timepoints)]
    temp_biased_barcodes <- get_biased_barcodes(temp, lineage_index = n_celltypes, lineage_comp_indices = 1:(n_celltypes-1),frac_cutoff=frac_cutoff,fold_cutoff=fold_cutoff)
    temp <- temp[temp_biased_barcodes, n_celltypes, drop = FALSE]
    temp <- reshape2::melt(as.matrix(temp))
    return(temp)
  })
  melted_data <- do.call(rbind, biased_barcodes_in_each_sample_only)
  colnames(melted_data) <- c("Barcode", "Sample", "Prop")
  melted_data$Sample <- factor(melted_data$Sample, levels = samples_for_plotting)
  thirtycolors <- c("#9843A1", "#CDAB44", "#A14854", "#8269DD", 
                    "#DD4469", "#A53F2A", "#6F8BDD", "#7D7078", "#C6D949", 
                    "#6497BF", "#80622E", "#CD91DC", "#D49980", "#D74FDB", 
                    "#E2442C", "#76D7D0", "#CBC6D0", "#5A8A34", "#584E8F", 
                    "#D77E33", "#C3D497", "#363550", "#D6489C", "#324026", 
                    "#6ADA44", "#532424", "#6AD88C", "#578B72", "#772B5D", 
                    "#D087A9","red","blue","green","gold","white","#9843A1", "#CDAB44", "#A14854", "#8269DD", 
                    "#DD4469", "#A53F2A", "#6F8BDD", "#7D7078", "#C6D949", 
                    "#6497BF", "#80622E", "#CD91DC", "#D49980", "#D74FDB", 
                    "#E2442C", "#76D7D0", "#CBC6D0", "#5A8A34", "#584E8F", 
                    "#D77E33", "#C3D497", "#363550", "#D6489C", "#324026", 
                    "#6ADA44", "#532424", "#6AD88C", "#578B72", "#772B5D", 
                    "#D087A9","red","blue","green","gold","white")
  set.seed(1)
  thirtycolors <- sample(thirtycolors)
  if(need_table){
    return(melted_data)  
  }
  
  for(i in 1:length(samples_for_plotting)){
    if(!samples_for_plotting[i]%in%unique(melted_data$Sample)){
      melted_data=rbind(melted_data,data.frame(Barcode="X",Sample=samples_for_plotting[i],Prop=0))
    }
  }

  ggplot2::ggplot(melted_data, ggplot2::aes(x = Sample, y = Prop, fill = Barcode))+
    ggplot2::scale_fill_manual(values = thirtycolors, guide = FALSE) +
    ggplot2::geom_bar(ggplot2::aes(y = Prop), stat = "identity", colour = "black") +
    ggplot2::scale_x_discrete("Months post-transplant",breaks=samples_for_plotting,labels=paste(timepoints,"m",sep="")) + 
    ggplot2::scale_y_continuous("Proportion",breaks = seq(from = 0, to = 1, by = 0.1), labels = function(x) {paste0(x * 100, "%")
    }, expand = c(0, 0)) + ggplot2::coord_cartesian(ylim = c(0,.75)) 
}

assign.bias=function(data,index=ncol(data)){
  other=data[,-index]
  selected=data[,index]
  
  max.other=apply(other,1,max)
  
  bias=selected/max.other
  
  bias=bias[!is.na(bias)]
  bias=as.data.frame(bias)
  
  barcodes.g10x=rownames(bias)[bias>=10]
  barcodes.5xg10x=rownames(bias)[bias<10&bias>=5]
  barcodes.2xg5x=rownames(bias)[bias<5&bias>=2]
  barcodes.1xg2x=rownames(bias)[bias<2&bias>=1]
  barcodes.1xg2xA=rownames(bias)[bias<1&bias>=1/2]
  barcodes.2xg5xA=rownames(bias)[bias<1/2&bias>=1/5]
  barcodes.5xg10xA=rownames(bias)[bias<1/5&bias>=1/10]
  barcodes.g10xA=rownames(bias)[bias<1/10]
  
  return(list(barcodes.g10x=barcodes.g10x,barcodes.5xg10x=barcodes.5xg10x,barcodes.2xg5x=barcodes.2xg5x,
              barcodes.1xg2x=barcodes.1xg2x,barcodes.1xg2xA=barcodes.1xg2xA,barcodes.2xg5xA=barcodes.2xg5xA,
              barcodes.5xg10xA=barcodes.5xg10xA,barcodes.g10xA=barcodes.g10xA))
}

bias.plot=function(data,index=ncol(data)){
  
  bias=assign.bias(data,index)
  
  barcodes.g10x=bias$barcodes.g10x
  barcodes.5xg10x=bias$barcodes.5xg10x
  barcodes.2xg5x=bias$barcodes.2xg5x
  barcodes.1xg2x=bias$barcodes.1xg2x
  
  df=NULL
  for(i in 1:4){
    barcodes=bias[[i]]
    values=data[barcodes,index]
    cat=c(">10x","5-10x","2-5x","1-2x")[i]
    if(length(barcodes)<1){
      barcodes=NA
      values=NA
    }
    temp.df=data.frame(barcodes=barcodes,values=values,cat=cat)
    df=rbind(df,temp.df)
  }
  
  df$color="black"
  df$color[df$values>.01]="red"
  
  
  if(sum(df$values,na.rm=TRUE)<1){
    print(ggplot(df,aes(x=cat,y=values,category=barcodes,fill=color))+
            geom_bar(color="black",stat="identity")+
            scale_fill_manual(values=c("gray48","red"))+
            scale_y_continuous("",limits=c(0,1),breaks=seq(0,1,by=.2),labels=c("0%","20%","40%","60%","80%","100%"))+
            scale_x_discrete("")+theme(legend.position="none")+theme_grey(base_size=22))
  }
  
  return(df)
  
}
assign.bias=function(data,index=ncol(data)){
  other=data[,-index]
  selected=data[,index]
  
  max.other=apply(other,1,max)
  
  bias=selected/max.other
  
  bias=bias[!is.na(bias)]
  bias=as.data.frame(bias)
  
  barcodes.g10x=rownames(bias)[bias>=10]
  barcodes.5xg10x=rownames(bias)[bias<10&bias>=5]
  barcodes.2xg5x=rownames(bias)[bias<5&bias>=2]
  barcodes.1xg2x=rownames(bias)[bias<2&bias>=1]
  barcodes.1xg2xA=rownames(bias)[bias<1&bias>=1/2]
  barcodes.2xg5xA=rownames(bias)[bias<1/2&bias>=1/5]
  barcodes.5xg10xA=rownames(bias)[bias<1/5&bias>=1/10]
  barcodes.g10xA=rownames(bias)[bias<1/10]
  
  return(list(barcodes.g10x=barcodes.g10x,barcodes.5xg10x=barcodes.5xg10x,barcodes.2xg5x=barcodes.2xg5x,
              barcodes.1xg2x=barcodes.1xg2x,barcodes.1xg2xA=barcodes.1xg2xA,barcodes.2xg5xA=barcodes.2xg5xA,
              barcodes.5xg10xA=barcodes.5xg10xA,barcodes.g10xA=barcodes.g10xA))
}

get_biased_barcodes <- function(your_data, lineage_index = ncol(your_data), lineage_comp_indices, frac_cutoff = 0.01, fold_cutoff = 10){
  your_data <- as.data.frame(prop.table(as.matrix(your_data), margin = 2))
  max_lineage_comp_indices <- apply(your_data[,lineage_comp_indices], 1, max)
  logical_lineage_comp_indices <- max_lineage_comp_indices*fold_cutoff < your_data[,lineage_index] & your_data[,lineage_index] > frac_cutoff
  return(rownames(your_data[logical_lineage_comp_indices,]))
}

cbc.all=function(cbc,day){
  ntime=nrow(cbc)
  df=NULL

  df$segs=cbc$WBC*cbc$SEGS/100
  df$rbc=cbc$RBC
  df$plt=cbc$PLT
  df$wbc=cbc$WBC
  df$grans=cbc$SEGS*cbc$WBC/100
  
  df$time=cbc$X.Daysinfus-day
  df=as.data.frame(df)
  df
  cbc.df=melt(df[,c("plt","rbc","segs")])
  cbc.df$time=rep(df$time,times=3)
  
  cbc.df$value=cbc.df$value*1000
  
  ggplot(cbc.df,aes(time,value,color=variable))+
    geom_point()+
    geom_line()+scale_x_continuous("Days post transplant",limits=c(-28,175),breaks=seq(-28,178,by=28))+
    scale_y_continuous("Number of cells per uL PB",trans = log10_trans(),breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x)))+ scale_colour_discrete(name="Cell type")+
    scale_color_manual(labels=c("PLT","RBC","SEGS"),values=c("red","green","blue"))
}


cbc.all2=function(cbc,day){
  ntime=nrow(cbc)
  df=NULL
  
  df$segs=cbc$WBC*cbc$SEGS/100
  df$rbc=cbc$RBC
  df$plt=cbc$PLT
  df$wbc=cbc$WBC
  df$grans=cbc$SEGS*cbc$WBC/100
  df$pbmnc=cbc$WBC*(cbc$LYMP+cbc$MONOS)/100
  df$monos=cbc$WBC*cbc$MONOS/100
  df$segs=cbc$WBC*cbc$SEGS/100
  df$lymp=cbc$WBC*cbc$LYMP/100
  df$cd3=cbc$CD3/cbc$P2*df$pbmnc
  df$cd20=cbc$CD20/cbc$P2*df$pbmnc
  
  
  df$time=cbc$X.Daysinfus-day
  df=as.data.frame(df)
  df
  cbc.df=melt(df[,c("cd20","monos")])
  cbc.df$time=rep(df$time,times=2)
  
  cbc.df$value=cbc.df$value*1000
  
  ggplot(cbc.df,aes(time,value,color=variable))+
    geom_point(size=3)+
    geom_line(size=2)+scale_x_continuous("Days post infection",limits=c(-28,175),breaks=seq(-28,178,by=28))+
    scale_y_continuous("Number of cells per uL PB",trans = log10_trans(),breaks = trans_breaks("log10", function(x) 10^x),limits=c(10^0,10^6),labels = trans_format("log10", math_format(10^.x)))+ scale_colour_discrete(name="Cell type")+
    scale_color_manual(labels=c("B cells","Monocytes"),values=c("red","green"))+theme_grey(base_size=22)
}

cbc.pre.zh33=function(cbc,day){
  ntime=nrow(cbc)
  df=NULL
  
  cbc$CD20=cbc$CD20/100
  cbc$CD3=cbc$CD3/100
  cbc$NK=cbc$NK/100
  
  df$pbmnc=cbc$WBC*(cbc$LYMP+cbc$MONOS)/100
  df$monos=cbc$WBC*cbc$MONOS/100
  df$segs=cbc$WBC*cbc$SEGS/100
  df$lymp=cbc$WBC*cbc$LYMP/100
  df$cd3=cbc$CD3/cbc$P2*df$pbmnc
  df$cd20=cbc$CD20/cbc$P2*df$pbmnc
  
  df$cd16=cbc$NK/cbc$P2*df$pbmnc
  #df$cd56=cbc$CD56/cbc$P2*df$pbmnc
  df$rbc=cbc$RBC
  df$plt=cbc$PLT
  df$wbc=cbc$WBC
  df$grans=cbc$SEGS*cbc$WBC/100
  
  df$time=cbc$X.Daysinfus-day
  df=as.data.frame(df)
  df
  cbc.df=melt(df[,c("pbmnc","cd20","cd3","cd16")])
  cbc.df$time=rep(df$time,times=4)
  
  cbc.df$value=cbc.df$value*1000
  
  # ggplot(cbc.df,aes(time,value,color=variable))+
  #   geom_point()+
  #   geom_line()+scale_x_continuous("Days post transplant",limits=c(0,175),breaks=seq(1,178,by=28))+
  #   scale_y_continuous("Number of cells",trans = log10_trans(),breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x)),limits=c(1,10^4))+ scale_colour_discrete(name="Cell type")+
  #   scale_color_manual(labels=c("PBMNC","B","T","NK"),values=c("black","gold","red","green"))
  return(cbc.df)
}

cbc.reinfection=function(cbc,day){
  ntime=nrow(cbc)
  df=NULL
  
  df$pbmnc=cbc$WBC*(cbc$LYMP+cbc$MONOS)/100
  df$monos=cbc$WBC*cbc$MONOS/100
  df$segs=cbc$WBC*cbc$SEGS/100
  df$lymp=cbc$WBC*cbc$LYMP/100
  df$cd3=cbc$CD3/cbc$P2*df$pbmnc
  df$cd8=cbc$CD8/cbc$P2*df$pbmnc
  
  df$cd4=cbc$CD4/cbc$P2*df$pbmnc
  df$cd20=cbc$CD20/cbc$P2*df$pbmnc
  
  df$nkg2a=cbc$NKG2A/cbc$P2*df$pbmnc
  df$cd16=cbc$CD16/cbc$P2*df$pbmnc
  df$cd56=cbc$CD56/cbc$P2*df$pbmnc
  df$rbc=cbc$RBC
  df$plt=cbc$PLT
  df$wbc=cbc$WBC
  df$grans=cbc$SEGS*cbc$WBC/100
  
  df$time=cbc$X.Daysinfus-day
  df=as.data.frame(df)
  
  cbc.df=melt(df[,c("cd4","cd8","cd16","cd56")])
  cbc.df$time=rep(df$time,times=4)
  
  cbc.df$value=cbc.df$value*1000
  print(ggplot(cbc.df,aes(time,value,color=variable))+
          geom_point(size=3)+
          geom_line(size=2)+scale_x_continuous("Days post re-rhCMV infection\n(Months post transplant)",limits=c(-3,30),breaks=c(seq(0,100,by=30))[1:2],labels=C(seq(0,100,by=30),paste("(",seq(19,22)," m)",sep=""))[1:2])+
          scale_y_continuous("Number of cells per uL PB",trans = log10_trans(),breaks = trans_breaks("log10", function(x) 10^x),limits=c(.01,10^6),labels = trans_format("log10", math_format(10^.x)))+ scale_colour_discrete(name="Cell type")+
          scale_color_manual("Cell Type",labels=c("T CD4+","T CD8+","NK CD56-CD16+","NK CD56+CD16-"),values=c("black","red","green","blue"))+theme_grey(base_size=22))
  return(cbc.df)
}


biased_barcode_summary=function(){
  
  df=data.frame(Monkey=NA,sample=NA,total=NA)
  
  jd76=biased_barcode_barplot(jd76.threshold[,jd76.all.samples.complete], length(jd76.time.complete), 5,jd76.time.complete, need_table = TRUE)
  for(i in unique(jd76$Sample)){
    df=rbind(df,c("JD76",i,sum(jd76$Prop[jd76$Sample==i])))
  }
  jm82=rbind(biased_barcode_barplot(jm82.threshold[,jm82.all.samples.complete], length(jm82.time.complete),5,jm82.time.complete, need_table = TRUE),
             c(NA,jm82.all.samples.complete[34],0))
  jm82=rbind(jm82,
             c(NA,jm82.all.samples.complete[35],0))
  for(i in unique(jm82$Sample)){
    df=rbind(df,c("JM82",i,sum(as.numeric(jm82$Prop)[jm82$Sample==i])))
  }
  zj31=biased_barcode_barplot(zj31.threshold[,zj31.all.samples.complete[-seq(5,25,by=5)]], length(zj31.time$cd56[-c(5,6)]), 5,zj31.time$cd56[-c(5,6)], need_table = TRUE)
  for(i in unique(zj31$Sample)){
    df=rbind(df,c("ZJ31",i,sum(zj31$Prop[zj31$Sample==i])))
  }
  zh33=biased_barcode_barplot(zh33.threshold[,zh33.all.samples.complete], length(zh33.time$nk), 4,zj31.time$nk, need_table = TRUE)
  for(i in unique(zh33$Sample)){
    df=rbind(df,c("ZH33",i,sum(zh33$Prop[zh33$Sample==i])))
  }
  zg66=biased_barcode_barplot(zg66.threshold[,zg66.all.samples], length(zg66.time$nk)-3, 4,zg66.time$nk[1:5], need_table = TRUE)
  for(i in unique(zg66$Sample)){
    df=rbind(df,c("ZG66",i,sum(zg66$Prop[zg66$Sample==i])))
  }
  zh19=biased_barcode_barplot(zh19.threshold[,zh19.all.samples], length(zh19.time$nk), 4,zh19.time$nk, need_table = TRUE)
  for(i in unique(zh19$Sample)){
    df=rbind(df,c("ZH19",i,sum(zh19$Prop[zh19$Sample==i])))
  }
  zk22=biased_barcode_barplot(zk22.threshold[,zk22.all.samples], length(zk22.time$nk[-4]), 5,zk22.time$nk[-4], need_table = TRUE)
  zk22=rbind(zk22,
             c(NA,zk22.all.samples[15],0))
  for(i in unique(zk22$Sample)){
    df=rbind(df,c("ZK22",i,sum(as.numeric(zk22$Prop)[zk22$Sample==i])))
  }
  
  df$sample=c(NA,jd76.time.complete,jm82.time.complete,zj31.time$cd56[-c(5,6)],zh33.time$nk[-2],zg66.time$nk[1:5],zh19.time$nk,zk22.time$nk[-4])
  df$total=as.numeric(df$total)
  
  df=df[-1,]
  ggplot(df,aes(sample,total,group=Monkey,color=Monkey))+
    geom_point(aes(shape=Monkey),size=4)+
    geom_line(size=1)+
    scale_color_manual(values=c("red","blue",rep("black",5)))+scale_y_continuous("Percent contribution of biased clones",breaks=seq(0,1,by=.1),labels=paste(seq(0,100,by=10),"%"))+
    scale_shape_manual(values=c(19,19,19,17,15,1,2,3))+scale_x_continuous("Months post transplant",limits=c(0,9.5),breaks=seq(1,10,by=2))+theme_gray(base_size = 22)
}

C=function(vector1,vector2){
  newvector=vector(length=length(vector1))
  for(i in 1:length(vector1)){
    newvector[i]=paste(vector1[i],"\n",vector2[i])
  }
  return(newvector)
}


autocorrelation=function(data,time){
  time=time[-1]
  x=cor(data,method="spearman")
  x=x[-1,]
  x=diag(x)
  x=data.frame(Autocorrelation=x,Months=time)
  
  print(x)
  
  print(ggplot(x,aes(x=Months,y=Autocorrelation))+
          geom_line(size=2)+
          geom_point(size=3)+scale_y_continuous(limits=c(0,1))+
          scale_x_continuous("Months post transplant",limits=c(9,20),breaks=seq(1,20,by=2))+
          theme_grey(base_size=22))
}

autocorrelation.cmvpos=function(Data,Time,Monkey){
  X=data.frame(Autocorrelation=NA,Months=NA,Monkey=NA)
  for(i in 1:length(Monkey)){
    time=Time[[i]]
    data=Data[[i]]
    #data=barcodetrackR::barcode_ggheatmap(data,printtable=TRUE,n_clones=500)
    time=time[-1]
    x=cor(data,method="spearman")
    x=x[-1,]
    x=diag(x)

    x=data.frame(Autocorrelation=x,Months=time,Monkey=Monkey[i])
    
    X=rbind(X,x)
  }

  
  print(ggplot(X[-1,],aes(x=Months,y=Autocorrelation,group=Monkey,color=Monkey))+
          geom_line(size=1)+
          geom_point(aes(shape=Monkey),size=4)+scale_y_continuous("Autocorrelation (Spearman)",limits=c(0,1))+
          scale_color_manual(values=c("red","blue",rep("black",5)))+
          scale_x_continuous("Months post transplant",limits=c(1,20),breaks=seq(1,20,by=2))+
          scale_shape_manual(values=c(19,19,19,17,15,1,2,3))+theme_grey(base_size=22))
}



barcode_ggheatmap_bar=function (your_data, names = colnames(your_data), n_clones = 10, 
                                your_title = "", grid = TRUE, label_size = 1, dendro = FALSE, 
                                cellnote_size = 4, printtable = FALSE, table_option = "percents", 
                                log_transform = TRUE, log_choice = exp(1), distance_method = "Euclidean", 
                                minkowski_power = 1, cellnote_option = "stars", hclust_linkage = "complete", 
                                row_order = "hierarchical", clusters = 0,
                                input="",COLORS=c("white","white","red","green")) 
{
  your_data_list <- list(raw_reads = your_data, prop_table = as.data.frame(prop.table(as.matrix(your_data), 
                                                                                      2)))
  your_data_list$prop_table[is.na(your_data)] <- 0
  if (any(colSums(your_data_list$prop_table) == 0)) {
    stop("One of your columns contained no data")
  }
  your_data_list$prop_table[your_data_list$prop_table == 0] <- NA
  your_data_list <- c(your_data_list, list(ranks = apply(-your_data_list$prop_table, 
                                                         2, rank, ties.method = "min", na.last = "keep")))
  top_clones_choices <- apply(your_data_list$ranks, 1, function(x) {
    any(x <= n_clones, na.rm = TRUE)
  })
  your_data_list <- lapply(your_data_list, function(x) {
    x[top_clones_choices, ]
  })
  your_data_list <- c(your_data_list, list(cellnote_matrix = your_data_list$ranks))
  your_data_list$cellnote_matrix[your_data_list$cellnote_matrix > 
                                   n_clones] <- NA
  your_data_list$cellnote_matrix[your_data_list$cellnote_matrix <= 
                                   n_clones] <- "*"
  your_data_list$prop_table[is.na(your_data_list$prop_table)] <- 0
  your_data_list <- c(your_data_list, list(logged_data = custom_log(your_data_list$prop_table, 
                                                                    log_choice, FALSE)))
  if (row_order == "hierarchical") {
    if (distance_method == "Minkowski") {
      hclustering <- hclust(proxy::dist((if (log_transform) 
        your_data_list$logged_data
        else your_data_list$prop_table), method = distance_method, 
        p = minkowski_power), method = hclust_linkage)
    }
    else {
      hclustering <- hclust(proxy::dist((if (log_transform) 
        your_data_list$logged_data
        else your_data_list$prop_table), method = distance_method), 
        method = hclust_linkage)
    }
    barcode_order <- rev(hclustering$order)
  }
  else if (row_order == "emergence") {
    barcode_order <- do.call(order, as.data.frame(your_data_list$prop_table))
  }
  else {
    stop("row_order must be one of \" hierarchical\" or \"emergence\"")
  }
  if (log_transform) {
    plotting_data <- reshape2::melt(as.matrix(your_data_list$logged_data))
    actual_scale <- c(log(100/4e+06, log_choice) - 1, log(100/4e+06, 
                                                          log_choice), log(0.001, log_choice), log(0.01, log_choice), 
                      log(0.1, log_choice), 0)
    your_scale <- scales::rescale(actual_scale, to = c(0, 
                                                       1))
    your_breaks <- c(actual_scale)
    your_labels <- c("Defined 0", "0.0025% ", "0.1%", "1%", 
                     "10%", "100%")
    your_limits <- c(log(100/4e+06, log_choice) - 1, 0)
    your_colors <- c("#4575B4", "#4575B4", "lightblue", "#fefeb9", 
                     "#D73027", "red4")
  }
  else {
    plotting_data <- reshape2::melt(as.matrix(your_data_list$prop_table))
    your_scale <- c(0, 100/4e+06, 0.005, 0.01, 0.1, 1)
    your_breaks <- c(0, 100/4e+06, 0.005, 0.01, 0.1, 1)
    your_labels <- c("0%", "0.0025%", "0.5%", "1%", "10%", 
                     "100%")
    your_limits <- c(0, 1)
    your_colors <- c("#4575B4", "#4575B4", "lightblue", "#fefeb9", 
                     "#D73027", "red4")
  }
  colnames(plotting_data) <- c("BARCODE", "SAMPLE", "Size")
  plotting_data=rbind(data.frame(BARCODE=row.names(input),SAMPLE="NA",Size=NA),plotting_data)
  plotting_data$color=c(as.vector(input[,1]),rep("",nrow(input)*(length(names)-1)))
  plotting_data$BARCODE <- factor(plotting_data$BARCODE, levels = rev(rownames(your_data_list$prop_table)[barcode_order]))
  if (cellnote_option %in% c("logs", "ranks", "stars", "percents", 
                             "reads")) {
    if (cellnote_option == "logs") {
      plotting_data$CELLNOTE <- round(reshape2::melt(as.matrix(your_data_list$logged_data))$value, 
                                      2)
      cellnote_vjust = 0.5
    }
    if (cellnote_option == "ranks") {
      plotting_data$CELLNOTE <- reshape2::melt(as.matrix(your_data_list$ranks))$value
      cellnote_vjust = 0.5
    }
    if (cellnote_option == "stars") {
      barcodes=length(unique(plotting_data$BARCODE))
      nsamples=length(unique(plotting_data$SAMPLE))
      plotting_data$CELLNOTE <- c(rep("*",barcodes),rep(NA,barcodes*(nsamples-1)))
      cellnote_vjust = 0.75
    }
    if (cellnote_option == "percents") {
      plotting_data$CELLNOTE <- paste0(round(reshape2::melt(as.matrix(your_data_list$prop_table))$value * 
                                               100, 2), "%")
      cellnote_vjust = 0.5
    }
    if (cellnote_option == "reads") {
      plotting_data$CELLNOTE <- reshape2::melt(as.matrix(your_data_list$raw_reads))$value
      cellnote_vjust = 0.5
    }
  }
  else {
    stop("cellnote_option must be one of c(\"logs\", \"ranks\", \"stars\", \"percents\", \"reads\")")
  }
  if (grid) 
    gridColor = "black"
  else gridColor = NA
  false_column_label <- unique(colnames(your_data_list$logged_data))[which(max(nchar(unique(colnames(your_data_list$logged_data)))) == 
                                                                             nchar(unique(colnames(your_data_list$logged_data))))[1]]
  if (printtable == TRUE) {
    switch(table_option, reads = return(your_data_list$raw_reads[barcode_order, 
                                                                 ]), percents = return(your_data_list$prop_table[barcode_order, 
                                                                                                                 ]), logs = return(your_data_list$logged_data[barcode_order, 
                                                                                                                                                              ]), ranks = return(your_data_list$ranks[barcode_order, 
                                                                                                                                                                                                      ]))
  }
  else {

    g1_heatmap <- ggplot2::ggplot(plotting_data, ggplot2::aes(x = SAMPLE, 
                                                              y = BARCODE)) + 
      ggplot2::geom_tile(ggplot2::aes(fill = Size),
                         color = gridColor) +
      ggplot2::scale_fill_gradientn(colors = your_colors, 
                                    values = your_scale, breaks = your_breaks, labels = your_labels, 
                                    limits = your_limits, expand = c(0, 0),na.value="white") + 
      ggplot2::scale_y_discrete(labels = NULL, 
                                breaks = NULL, expand = c(0, 0)) + ggplot2::scale_x_discrete(expand = c(0, 
                                                                                                        0), 
                                                                                             labels = names) + 
      ggplot2::ylab(NULL) + ggplot2::xlab(NULL) + 
      ggplot2::ggtitle(paste0("\n", your_title, "\n")) + 
      ggplot2::geom_text(ggplot2::aes(label = CELLNOTE,colour=color), 
                         vjust = cellnote_vjust, size = cellnote_size) + 
      scale_colour_manual(values=COLORS)+
      ggplot2::theme(plot.title = ggplot2::element_text(size = 20), 
                     axis.text.x = ggplot2::element_text(angle = 90, 
                                                         hjust = 1, vjust = 0.5, size = label_size), 
                     legend.text = ggplot2::element_text(size = 15, 
                                                         face = "bold"), legend.title = ggplot2::element_text(size = 15), 
                     axis.ticks = ggplot2::element_blank(), legend.key.height = ggplot2::unit(5, 
                                                                                              "line"))
    if (row_order == "emergence") {
      g1_heatmap
    }
    else {
      dendro_data <- ggdendro::dendro_data(hclustering, 
                                           type = "rectangle")
      g2_dendrogram <- ggplot2::ggplot(ggdendro::segment(dendro_data)) + 
        ggplot2::geom_segment(ggplot2::aes(x = x, y = y, 
                                           xend = xend, yend = yend)) + ggplot2::scale_x_discrete(expand = c(0.5/nrow(your_data_list$logged_data), 
                                                                                                             0.01)) + ggplot2::scale_y_reverse(expand = c(0.01, 
                                                                                                                                                          0), labels = false_column_label, breaks = if (row_order == 
                                                                                                                                                                                                        "emergence") 
                                                                                                                                                            NULL
                                                                                                                                               else mean(hclustering$height)) + ggplot2::coord_flip() + 
        ggplot2::ylab(NULL) + ggplot2::xlab(NULL) + ggplot2::ggtitle("\n \n") + 
        ggplot2::theme(plot.title = ggplot2::element_text(size = 20), 
                       axis.text.x = ggplot2::element_text(colour = "white", 
                                                           angle = 90, hjust = 1, vjust = 0.5, size = label_size), 
                       panel.background = ggplot2::element_rect(fill = "white", 
                                                                colour = "white"), axis.ticks = ggplot2::element_blank())
      if (dendro) {
        if (clusters > 0) {
          custom_colors <- c("#89C5DA", "#DA5724", "#74D944", 
                             "#CE50CA", "#3F4921", "#C0717C", "#CBD588")
          clusters_df <- data.frame(CLUSTERS = cutree(hclustering, 
                                                      clusters), ASSIGNMENT = factor(hclustering$labels, 
                                                                                     levels = hclustering$labels[(hclustering$order)]))
          g3_clusters <- ggplot2::ggplot(clusters_df, 
                                         ggplot2::aes(x = 1, y = ASSIGNMENT, fill = factor(CLUSTERS))) + 
            ggplot2::geom_tile() + ggplot2::ggtitle("\n \n") + 
            ggplot2::scale_fill_manual(values = custom_colors[1:clusters]) + 
            ggplot2::scale_x_continuous(expand = c(0, 
                                                   0), labels = false_column_label, breaks = 1) + 
            ggplot2::theme(plot.title = ggplot2::element_text(size = 20), 
                           axis.title = ggplot2::element_blank(), 
                           axis.ticks = ggplot2::element_blank(), 
                           axis.text.y = ggplot2::element_blank(), 
                           axis.text.x = ggplot2::element_text(colour = "white", 
                                                               angle = 90, hjust = 1, vjust = 0.5, size = label_size), 
                           legend.position = "none")
          gridExtra::grid.arrange(g2_dendrogram, g3_clusters, 
                                  g1_heatmap, ncol = 3, widths = c(1, 0.2, 
                                                                   4))
        }
        else {
          gridExtra::grid.arrange(g2_dendrogram, g1_heatmap, 
                                  ncol = 2, widths = c(1, 4))
        }
      }
      else {
        g1_heatmap
      }
    }
  }
}


barcode_ggheatmapN_bar=function (your_data,your_facs_data, names = colnames(your_data), n_clones = 10, 
                                your_title = "", grid = TRUE, label_size = 1, dendro = FALSE, 
                                cellnote_size = 4, printtable = FALSE, table_option = "percents", 
                                log_transform = TRUE, log_choice = exp(1), distance_method = "Euclidean", 
                                minkowski_power = 1, cellnote_option = "stars", hclust_linkage = "complete", 
                                row_order = "hierarchical", clusters = 0,
                                input="",COLORS=c("white","white","red","black")) 
{
  
  your_data_list <- list(raw_reads = your_data, prop_table = as.data.frame(prop.table(as.matrix(your_data), 
                                                                                      2)))
  your_data_list$prop_table=t(t(your_data_list$prop_table)*your_facs_data)
  your_data_list$prop_table[is.na(your_data)] <- 0
  if (any(colSums(your_data_list$prop_table) == 0)) {
    stop("One of your columns contained no data")
  }
  your_data_list$prop_table[your_data_list$prop_table == 0] <- NA
  your_data_list <- c(your_data_list, list(ranks = apply(-your_data_list$prop_table, 
                                                         2, rank, ties.method = "min", na.last = "keep")))
  top_clones_choices <- apply(your_data_list$ranks, 1, function(x) {
    any(x <= n_clones, na.rm = TRUE)
  })
  your_data_list <- lapply(your_data_list, function(x) {
    x[top_clones_choices, ]
  })
  your_data_list <- c(your_data_list, list(cellnote_matrix = your_data_list$ranks))
  your_data_list$cellnote_matrix[your_data_list$cellnote_matrix > 
                                   n_clones] <- NA
  your_data_list$cellnote_matrix[your_data_list$cellnote_matrix <= 
                                   n_clones] <- "*"
  your_data_list$prop_table[is.na(your_data_list$prop_table)] <- 0
  your_data_list <- c(your_data_list, list(logged_data = custom_log(your_data_list$prop_table, 
                                                                    log_choice, FALSE)))
  if (row_order == "hierarchical") {
    if (distance_method == "Minkowski") {
      hclustering <- hclust(proxy::dist((if (log_transform) 
        your_data_list$logged_data
        else your_data_list$prop_table), method = distance_method, 
        p = minkowski_power), method = hclust_linkage)
    }
    else {
      hclustering <- hclust(proxy::dist((if (log_transform) 
        your_data_list$logged_data
        else your_data_list$prop_table), method = distance_method), 
        method = hclust_linkage)
    }
    barcode_order <- rev(hclustering$order)
  }
  else if (row_order == "emergence") {
    barcode_order <- do.call(order, as.data.frame(your_data_list$prop_table))
  }
  else {
    stop("row_order must be one of \" hierarchical\" or \"emergence\"")
  }
  if (log_transform) {
    plotting_data <- reshape2::melt(as.matrix(your_data_list$logged_data))
    actual_scale <- c(log(100/4e+06, log_choice) - 1, log(100/4e+06, 
                                                          log_choice), log(0.001, log_choice), log(0.01, log_choice), 
                      log(0.1, log_choice), 0)
    your_scale <- scales::rescale(actual_scale, to = c(0, 
                                                       1))
    your_breaks <- c(actual_scale)
    your_labels <- c("Defined 0", "0.0025% ", "0.1%", "1%", 
                     "10%", "100%")
    your_limits <- c(log(100/4e+06, log_choice) - 1, 0)
    your_colors <- c("#4575B4", "#4575B4", "lightblue", "#fefeb9", 
                     "#D73027", "red4")
  }
  else {
    plotting_data <- reshape2::melt(as.matrix(your_data_list$prop_table))
    your_scale <- c(0, 100/4e+06, 0.005, 0.01, 0.1, 1)
    your_breaks <- c(0, 100/4e+06, 0.005, 0.01, 0.1, 1)
    your_labels <- c("0%", "0.0025%", "0.5%", "1%", "10%", 
                     "100%")
    your_limits <- c(0, 1)
    your_colors <- c("#4575B4", "#4575B4", "lightblue", "#fefeb9", 
                     "#D73027", "red4")
  }
  colnames(plotting_data) <- c("BARCODE", "SAMPLE", "Size")
  plotting_data=rbind(data.frame(BARCODE=row.names(input),SAMPLE="NA",Size=NA),plotting_data)
  plotting_data$color=c(as.vector(input[,1]),rep("",nrow(input)*(length(names)-1)))
  plotting_data$BARCODE <- factor(plotting_data$BARCODE, levels = rev(rownames(your_data_list$prop_table)[barcode_order]))
  if (cellnote_option %in% c("logs", "ranks", "stars", "percents", 
                             "reads")) {
    if (cellnote_option == "logs") {
      plotting_data$CELLNOTE <- round(reshape2::melt(as.matrix(your_data_list$logged_data))$value, 
                                      2)
      cellnote_vjust = 0.5
    }
    if (cellnote_option == "ranks") {
      plotting_data$CELLNOTE <- reshape2::melt(as.matrix(your_data_list$ranks))$value
      cellnote_vjust = 0.5
    }
    if (cellnote_option == "stars") {
      barcodes=length(unique(plotting_data$BARCODE))
      nsamples=length(unique(plotting_data$SAMPLE))
      plotting_data$CELLNOTE <- c(rep("*",barcodes),rep(NA,barcodes*(nsamples-1)))
      cellnote_vjust = 0.75
    }
    if (cellnote_option == "percents") {
      plotting_data$CELLNOTE <- paste0(round(reshape2::melt(as.matrix(your_data_list$prop_table))$value * 
                                               100, 2), "%")
      cellnote_vjust = 0.5
    }
    if (cellnote_option == "reads") {
      plotting_data$CELLNOTE <- reshape2::melt(as.matrix(your_data_list$raw_reads))$value
      cellnote_vjust = 0.5
    }
  }
  else {
    stop("cellnote_option must be one of c(\"logs\", \"ranks\", \"stars\", \"percents\", \"reads\")")
  }
  if (grid) 
    gridColor = "black"
  else gridColor = NA
  false_column_label <- unique(colnames(your_data_list$logged_data))[which(max(nchar(unique(colnames(your_data_list$logged_data)))) == 
                                                                             nchar(unique(colnames(your_data_list$logged_data))))[1]]
  if (printtable == TRUE) {
    switch(table_option, reads = return(your_data_list$raw_reads[barcode_order, 
                                                                 ]), percents = return(your_data_list$prop_table[barcode_order, 
                                                                                                                 ]), logs = return(your_data_list$logged_data[barcode_order, 
                                                                                                                                                              ]), ranks = return(your_data_list$ranks[barcode_order, 
                                                                                                                                                                                                      ]))
  }
  else {
    
    g1_heatmap <- ggplot2::ggplot(plotting_data, ggplot2::aes(x = SAMPLE, 
                                                              y = BARCODE)) + 
      ggplot2::geom_tile(ggplot2::aes(fill = Size),
                         color = NA) +
      ggplot2::scale_fill_gradient(low="white",high="purple4",na.value="white")+
      ggplot2::scale_y_discrete(labels = NULL, 
                                breaks = NULL, expand = c(0, 0)) + ggplot2::scale_x_discrete(expand = c(0, 
                                                                                                        0), 
                                                                                             labels = names) + 
      ggplot2::ylab(NULL) + ggplot2::xlab(NULL) + 
      ggplot2::ggtitle(paste0("\n", your_title, "\n")) + 
      ggplot2::geom_text(ggplot2::aes(label = CELLNOTE,colour=color), 
                         vjust = cellnote_vjust, size = cellnote_size) + 
      scale_colour_manual(values=COLORS)+
      ggplot2::theme(plot.title = ggplot2::element_text(size = 20), 
                     axis.text.x = ggplot2::element_text(angle = 90, 
                                                         hjust = 1, vjust = 0.5, size = label_size), 
                     legend.text = ggplot2::element_text(size = 15, 
                                                         face = "bold"), legend.title = ggplot2::element_text(size = 15), 
                     axis.ticks = ggplot2::element_blank(), legend.key.height = ggplot2::unit(5, 
                                                                                              "line"))
    if (row_order == "emergence") {
      g1_heatmap
    }
    else {
      dendro_data <- ggdendro::dendro_data(hclustering, 
                                           type = "rectangle")
      g2_dendrogram <- ggplot2::ggplot(ggdendro::segment(dendro_data)) + 
        ggplot2::geom_segment(ggplot2::aes(x = x, y = y, 
                                           xend = xend, yend = yend)) + ggplot2::scale_x_discrete(expand = c(0.5/nrow(your_data_list$logged_data), 
                                                                                                             0.01)) + ggplot2::scale_y_reverse(expand = c(0.01, 
                                                                                                                                                          0), labels = false_column_label, breaks = if (row_order == 
                                                                                                                                                                                                        "emergence") 
                                                                                                                                                            NULL
                                                                                                                                               else mean(hclustering$height)) + ggplot2::coord_flip() + 
        ggplot2::ylab(NULL) + ggplot2::xlab(NULL) + ggplot2::ggtitle("\n \n") + 
        ggplot2::theme(plot.title = ggplot2::element_text(size = 20), 
                       axis.text.x = ggplot2::element_text(colour = "white", 
                                                           angle = 90, hjust = 1, vjust = 0.5, size = label_size), 
                       panel.background = ggplot2::element_rect(fill = "white", 
                                                                colour = "white"), axis.ticks = ggplot2::element_blank())
      if (dendro) {
        if (clusters > 0) {
          custom_colors <- c("#89C5DA", "#DA5724", "#74D944", 
                             "#CE50CA", "#3F4921", "#C0717C", "#CBD588")
          clusters_df <- data.frame(CLUSTERS = cutree(hclustering, 
                                                      clusters), ASSIGNMENT = factor(hclustering$labels, 
                                                                                     levels = hclustering$labels[(hclustering$order)]))
          g3_clusters <- ggplot2::ggplot(clusters_df, 
                                         ggplot2::aes(x = 1, y = ASSIGNMENT, fill = factor(CLUSTERS))) + 
            ggplot2::geom_tile() + ggplot2::ggtitle("\n \n") + 
            ggplot2::scale_fill_manual(values = custom_colors[1:clusters]) + 
            ggplot2::scale_x_continuous(expand = c(0, 
                                                   0), labels = false_column_label, breaks = 1) + 
            ggplot2::theme(plot.title = ggplot2::element_text(size = 20), 
                           axis.title = ggplot2::element_blank(), 
                           axis.ticks = ggplot2::element_blank(), 
                           axis.text.y = ggplot2::element_blank(), 
                           axis.text.x = ggplot2::element_text(colour = "white", 
                                                               angle = 90, hjust = 1, vjust = 0.5, size = label_size), 
                           legend.position = "none")
          gridExtra::grid.arrange(g2_dendrogram, g3_clusters, 
                                  g1_heatmap, ncol = 3, widths = c(1, 0.2, 
                                                                   4))
        }
        else {
          gridExtra::grid.arrange(g2_dendrogram, g1_heatmap, 
                                  ncol = 2, widths = c(1, 4))
        }
      }
      else {
        g1_heatmap
      }
    }
  }
}

diversity=function(data,time){
  diversity=NULL
  for(i in colnames(data)){
    temp=data[,i]
    temp.df=data.frame(Month=i,Count=temp,Barcode=rownames(data))
    temp_diversity=vegan::diversity(data[,i])
    
    diversity=c(diversity,as.numeric(temp_diversity))
  }
  diversity=data.frame(Diversity=diversity,Cell=c(rep("Grans",length(time$grans)),rep("T",length(time$t)),rep("B",length(time$b)),rep("NK",length(time$nk))),time=(melt(time)$value))
  
  ggplot(diversity,aes(x=time,y=Diversity,color=Cell))+
    geom_line(size=2)+
    geom_point(size=3)+theme_grey(base_size=22)+scale_y_continuous(limits=c(4,8.5))+
    scale_x_continuous("Months post transplant",breaks=seq(1,9.5,by=2),limits=c(1,10))
}

diversity2=function(data,time){
  diversity=NULL
  for(i in colnames(data)){
    temp=data[,i]
    temp.df=data.frame(Month=i,Count=temp,Barcode=rownames(data))
    temp_diversity=vegan::diversity(data[,i])
    
    diversity=c(diversity,as.numeric(temp_diversity))
  }
  diversity=data.frame(Diversity=diversity,Cell=c(rep("Grans",length(time$grans)),rep("T",length(time$t)),rep("B",length(time$b)),rep("NK CD56-CD16+",length(time$nk)),rep("NK CD56+CD16-",length(time$cd56))),time=(melt(time)$value))
  
  ggplot(diversity,aes(x=time,y=Diversity,color=Cell))+
    geom_line(size=2)+
    geom_point(size=3)+theme_grey(base_size=22)+scale_y_continuous(limits=c(4,8.5))+
    scale_x_continuous("Months post transplant",breaks=seq(1,9.5,by=2),limits=c(1,10))
}

