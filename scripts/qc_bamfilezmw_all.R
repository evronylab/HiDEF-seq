#!/usr/bin/env Rscript
# Usage: qc_bamfilezmw_all.R [configuration.yaml]
#  - Uses [configuration.yaml] file per documentation
#  - Outputs QC plots of bamfilezmw_all
#
# Prerequisite R packages: ggplot2, plyr, stringr, configr, qs
#
######################
### Load required libraries
######################
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(configr))
suppressPackageStartupMessages(library(qs))

######################
### Load configuration
######################
args <- commandArgs(trailingOnly = TRUE)
if (length(args)!=1) {
  stop("Usage: Rscript qc_bamfilezmw_all.R [configuration.yaml]\n", call.=FALSE)
}
cat("#### Loading configuration...")
yaml.config <- suppressWarnings(read.config(args[1]))

bamfilezmw_all_filename <- yaml.config$bamfilezmw_all_filename
outputpath <- paste0(yaml.config$analysisoutput_path,"/")

if(!dir.exists(outputpath)){dir.create(outputpath,recursive=TRUE)}

cat("DONE\n")

cat("#### Loading bamfilezmw_all...")
bamfilezmw_all <- qread(bamfilezmw_all_filename)
cat("DONE\n")

#Histograms of rq (alignment quality), ec (effective coverage), np (number of passes). Plot for each sample, and plot fwd, rev, and avg of ccs reads. Only need to plot for one reference genome because these parameters don't depend on the reference genome. Also add median value to plot.
for(i in names(bamfilezmw_all)){
  for(j in names(bamfilezmw_all[[i]])) {
  	if(class(bamfilezmw_all[[i]][[j]]) != "list"){
  		#Skip this sample
  		cat(" No reads or mutations in sample",paste0(i,".",j),"in selected chromosomes. Proceeding to next sample...\n")
  		next
  	}
  		
    for(k in c("rq","ec","np")){
      histvalues.fwd <- bamfilezmw_all[[i]][[j]]$tag.fwd[[k]]
      histvalues.rev <- bamfilezmw_all[[i]][[j]]$tag.rev[[k]]
      histvalues.avg <- apply(cbind(histvalues.fwd,histvalues.rev),1,mean)
      df <- data.frame(label=c(rep("fwd",length(histvalues.fwd)),rep("rev",length(histvalues.rev)),rep("fwdrevavg",length(histvalues.avg))),values=c(histvalues.fwd,histvalues.rev,histvalues.avg))
      df$label <- factor(df$label, levels = unique(df$label))
      ggplot(df, aes(x=values,color=label)) + geom_density(adjust=1/2) + labs(x=k,subtitle=paste("Median of fwd/rev avg values:",median(histvalues.avg)))
      ggsave(paste0(outputpath,i,".",j,".",k,".hist.pdf"))
    }
  }
}

#Histogram of qwidth, isize and mapq (fwd, rev, and average of fwd and rev ccs). Plot for each sample and for each reference genome. Also add median value of isize to plot.
for(i in names(bamfilezmw_all)){
  for(j in names(bamfilezmw_all[[i]])) {
    if(class(bamfilezmw_all[[i]][[j]]) != "list"){
  		#Skip this sample
  		cat(" No reads or mutations in sample",paste0(i,".",j),"in selected chromosomes. Proceeding to next sample...\n")
  		next
    }
    histvalues.fwd <- bamfilezmw_all[[i]][[j]]$qwidth.fwd
    histvalues.rev <- bamfilezmw_all[[i]][[j]]$qwidth.rev
    histvalues.avg <- apply(cbind(histvalues.fwd,histvalues.rev),1,mean)
    df <- rbind(data.frame(label=c(rep("fwd",length(histvalues.fwd)),rep("rev",length(histvalues.rev)),rep("fwdrevavg",length(histvalues.avg))),values=c(histvalues.fwd,histvalues.rev,histvalues.avg)))
    df$label <- factor(df$label, levels = unique(df$label))
    ggplot(df, aes(x=values,color=label)) + geom_density(adjust=1/2) + labs(x="qwidth",subtitle=paste("Median of fwd/rev avg values:",median(histvalues.avg)))
    ggsave(paste0(outputpath,i,".",j,".qwidth.hist.pdf"))
    
    histvalues.fwd <- bamfilezmw_all[[i]][[j]]$isize.fwd
    histvalues.rev <- bamfilezmw_all[[i]][[j]]$isize.rev
    histvalues.avg <- apply(cbind(histvalues.fwd,histvalues.rev),1,mean)
    df <- rbind(data.frame(label=c(rep("fwd",length(histvalues.fwd)),rep("rev",length(histvalues.rev)),rep("fwdrevavg",length(histvalues.avg))),values=c(histvalues.fwd,histvalues.rev,histvalues.avg)))
    df$label <- factor(df$label, levels = unique(df$label))
    ggplot(df, aes(x=values,color=label)) + geom_density(adjust=1/2) + labs(x="isize",subtitle=paste("Median of fwd/rev avg values:",median(histvalues.avg)))
    ggsave(paste0(outputpath,i,".",j,".isize.hist.pdf"))
  
    histvalues.fwd <- bamfilezmw_all[[i]][[j]]$mapq.fwd
    histvalues.rev <- bamfilezmw_all[[i]][[j]]$mapq.rev
    histvalues.avg <- apply(cbind(histvalues.fwd,histvalues.rev),1,mean)
    df <- rbind(data.frame(label=c(rep("fwd",length(histvalues.fwd)),rep("rev",length(histvalues.rev)),rep("fwdrevavg",length(histvalues.avg))),values=c(histvalues.fwd,histvalues.rev,histvalues.avg)))
    df$label <- factor(df$label, levels = unique(df$label))
    df.new <- ddply(df,.(label),summarise,
          prop=as.numeric(prop.table(table(values))),
          values=as.numeric(names(table(values))))
    ggplot(df.new,aes(x=values,y=prop,fill=label)) + geom_col(position="dodge") + labs(x="mapq") + ylim(0,1)
    #+ xlim(50,61)
    #ggplot(df, aes(x=values,color=label)) + geom_density(adjust=1/2) + labs(x="mapq")
    ggsave(paste0(outputpath,i,".",j,".mapq.hist.pdf")) 
  }
}

#Scatter plots between rq, ec, np, isize, mapq (all combinations, for each take avg of fwd and rev). Plot for each reference genome. Note: to plot log scale of counts add aes(fill=log(..count..)) to geom() function.
tagsforscatterplot <- c("rq","ec","np","qwidth","isize","mapq")
for(i in names(bamfilezmw_all)){
  for(r in names(bamfilezmw_all[[i]])){
  	if(class(bamfilezmw_all[[i]][[r]]) != "list"){
  		#Skip this sample
  		cat(" No reads or mutations in sample",paste0(i,".",r),"in selected chromosomes. Proceeding to next sample...\n")
  		next
  	}
  	
    for(j in tagsforscatterplot){
      for(k in setdiff(tagsforscatterplot,j)){
        
        if(j == "isize"){
          avgj=apply(cbind(bamfilezmw_all[[i]][[r]]$isize.fwd,bamfilezmw_all[[i]][[r]]$isize.rev),1,mean)
        }else if (j=="qwidth"){
          avgj=apply(cbind(bamfilezmw_all[[i]][[r]]$qwidth.fwd,bamfilezmw_all[[i]][[r]]$qwidth.rev),1,mean)
        }else if (j=="mapq"){
          avgj=apply(cbind(bamfilezmw_all[[i]][[r]]$mapq.fwd,bamfilezmw_all[[i]][[r]]$mapq.rev),1,mean)
        }else{
          avgj=apply(cbind(bamfilezmw_all[[i]][[r]]$tag.fwd[[j]],bamfilezmw_all[[i]][[r]]$tag.rev[[j]]),1,mean)
        }
        
        if(k == "isize"){
          avgk=apply(cbind(bamfilezmw_all[[i]][[r]]$isize.fwd,bamfilezmw_all[[i]][[r]]$isize.rev),1,mean)
        }else if (k=="qwidth"){
          avgk=apply(cbind(bamfilezmw_all[[i]][[r]]$qwidth.fwd,bamfilezmw_all[[i]][[r]]$qwidth.rev),1,mean)
        }else if (k=="mapq"){
          avgk=apply(cbind(bamfilezmw_all[[i]][[r]]$mapq.fwd,bamfilezmw_all[[i]][[r]]$mapq.rev),1,mean)
        }else{
          avgk=apply(cbind(bamfilezmw_all[[i]][[r]]$tag.fwd[[k]],bamfilezmw_all[[i]][[r]]$tag.rev[[k]]),1,mean)
        }
        
        df <- data.frame(avgjdata=avgj,avgkdata=avgk)
        ggplot(df,aes(x=avgjdata,y=avgkdata)) + geom_bin2d(bins=250) + labs(x=j,y=k)
        ggsave(paste0(outputpath,i,".",r,".",j,".",k,".scatter.pdf"))
        
      }
    }
  }
}

#Scatter plots of fwd vs rev for rq, ec, np, qwidth, mapq, isize.
for(i in names(bamfilezmw_all)){
  for(j in names(bamfilezmw_all[[i]])){
  	if(class(bamfilezmw_all[[i]][[j]]) != "list"){
  		#Skip this sample
  		cat(" No reads or mutations in sample",paste0(i,".",j),"in selected chromosomes. Proceeding to next sample...\n")
  		next
  	}
  	
    for(k in c("rq","ec","np")){
    df <- data.frame(fwd=bamfilezmw_all[[i]][[j]]$tag.fwd[[k]],rev=bamfilezmw_all[[i]][[j]]$tag.rev[[k]])
    ggplot(df,aes(x=fwd,y=rev)) + geom_bin2d(bins=250) + ggtitle(paste(i,j,k,"fwd vs rev scatter plot"))
    ggsave(paste0(outputpath,i,".",j,".",k,".","fwdvsrev.scatter.pdf"))
    }
  }
}

for(i in names(bamfilezmw_all)){
  for(j in names(bamfilezmw_all[[i]])){
  	if(class(bamfilezmw_all[[i]][[j]]) != "list"){
  		#Skip this sample
  		cat(" No reads or mutations in sample",paste0(i,".",j),"in selected chromosomes. Proceeding to next sample...\n")
  		next
  	}
  	
    for(k in c("isize","qwidth","mapq")){
    if(k=="isize" | k=="qwidth"){
      df <- data.frame(fwd=unlist(bamfilezmw_all[[i]][[j]][[paste0(k,".fwd")]]),rev=unlist(bamfilezmw_all[[i]][[j]][[paste0(k,".rev")]]))
      ggplot(df,aes(x=fwd,y=rev)) + geom_bin2d(bins=250) + ggtitle(paste(i,j,k,"fwd vs rev scatter plot"))
    }else if(k=="mapq"){
      df <- data.frame(fwd=bamfilezmw_all[[i]][[j]][[paste0(k,".fwd")]],rev=bamfilezmw_all[[i]][[j]][[paste0(k,".rev")]])
      ggplot(df,aes(x=fwd,y=rev)) + geom_bin2d(bins=60,color="black") + scale_fill_gradient(low = "white", high = "black") + theme_light() + ggtitle(paste(i,j,k,"fwd vs rev scatter plot"))
    }
    ggsave(paste0(outputpath,i,".",j,".",k,".","fwdvsrev.scatter.pdf"))
    }
  }
}

#Histogram of vcf Depth, VAF, GQ (fwd, rev, and ccs mutations). Plot for each sample, reference genome, and VCF file.
for(i in names(bamfilezmw_all)){
  for(j in names(bamfilezmw_all[[i]])) {
  	if(class(bamfilezmw_all[[i]][[j]]) != "list"){
  		#Skip this sample
  		cat(" No reads or mutations in sample",paste0(i,".",j),"in selected chromosomes. Proceeding to next sample...\n")
  		next
  	}
  	
    for(q in names(bamfilezmw_all[[i]][[j]]$vcfSNV.fwd))
      for(k in c("Depth","VAF","GQ")){
        histvalues.fwd <- unlist(lapply(bamfilezmw_all[[i]][[j]]$vcfSNV.fwd[[q]][[k]],function(x){x[x>=0 & !is.na(x)]}))
        histvalues.rev <- unlist(lapply(bamfilezmw_all[[i]][[j]]$vcfSNV.rev[[q]][[k]],function(x){x[x>=0 & !is.na(x)]}))
        histvalues.zmw <- unlist(lapply(bamfilezmw_all[[i]][[j]]$zmw.vcfSNV[[q]][[k]],function(x){x[x>=0 & !is.na(x)]}))
        
        df <- rbind(data.frame(label=c(rep("fwd",length(histvalues.fwd)),rep("rev",length(histvalues.rev)),rep("shared",length(histvalues.zmw))),values=c(histvalues.fwd,histvalues.rev,histvalues.zmw)))
        df$label <- factor(df$label, levels = unique(df$label))
        df$values <- as.numeric(df$values)
        
        if(k=="Depth"){
          ggplot(df, aes(x=values,color=label)) + geom_density(adjust=1/2) + labs(x=k) + scale_x_continuous(trans="log10")
        }else{
          ggplot(df, aes(x=values,color=label)) + geom_density(adjust=1/2) + labs(x=k)
        }
        ggsave(paste0(outputpath,i,".",j,".",basename(q),".",k,".hist.pdf"))
      }
  }
}
cat("\nDONE\n")
