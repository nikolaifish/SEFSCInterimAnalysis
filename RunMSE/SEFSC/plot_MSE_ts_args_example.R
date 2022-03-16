OMNames= OMNames
scenNames=scenNames_user
myvarName="Catch"
myvarDenomExpr = NULL
ylabel = myvarName
FUN = mean
#sims = NULL # which simulation runs should be plotted?
sims = 1:3
MPNamesSub = "SCA_1"#MPNamesSub_user
MPNamesSub_legtext="Annual assessment"#MPNamesSub_user_legtext
rowLabs=scenNames_key[scenNames_user]
MPNamesRef = NULL # Name of reference MP
FUNRef = function(myvarsubStat,MPNamesRef){ # Function for comparing MPs with reference MP
  (myvarsubStat-myvarsubStat[,MPNamesRef])/myvarsubStat[,MPNamesRef]}
legtext_o = NULL # order for plotting legend text
MSEoutNames = paste0(rep(gsub("OM","MSE",OMNames),each=length(scenNames)),
                     "_",
                     rep(scenNames,length(OMNames))
)
colLabs = str_replace_all(gsub("OM_","",OMNames),"(?<=[a-z])(?=[A-Z])"," ")
plot_type= "sims"#"huynh" # "huynh" plots time series by MP like in the Huynh et al paper
# "boot" plots times series of bootstrapped data  with bamExtras::plot_boot_vec
# "phase" plots phase plots of F/Fmsy versus SSB/SSBmsy
# "tradeoff" plots tradeoff plots by MP
# "addind" plot additional indices from AddInd
# "sims" plot time series for individual sim runs
proyear_phase=NULL
PMx = Yield # Performance metric function x (x-axis) for tradeoff plots)
PMy = P75 # Performance metric function y (y-axis) for tradeoff plots)
mar=c(1,1,1,1)
oma=c(6,4,3,3)
mgp=c(1.5,0.2,0)
cols = setNames(ROGBP(length(MPNamesSub)),MPNamesSub)
ltys = setNames(rep(1,length(MPNamesSub)),MPNamesSub)
pchs = setNames(rep(16,length(MPNamesSub)),MPNamesSub)
lwds = setNames(rep(3,length(MPNamesSub)),MPNamesSub)
types = setNames(rep("l",length(MPNamesSub)),MPNamesSub)
xlabel = "Management year"
xlim = NULL
ylim = NULL
assessmentYears = seq(10,40,by=10)
x.at = NULL
y.at = NULL
# leg.x=NULL
# leg.y =NULL
legend_panel=NULL
legend_cex= 1
legend_x="bottomleft"
legend_inset=c(0,-0.6)
label_axes="all" # "all" labels all panels, "leftbottom" labels y axes on the left and x axes on the bottom panels
hline=1 # values to pass into abline(h=hline)

  if(is.null(legtext_o)){
    legtext_o <- seq_along(MPNamesSub)
  }
  
  xlim_user <- xlim
  ylim_user <- ylim
  x.at_user <- x.at
  y.at_user <- y.at
  
  colLabsAll = as.character(matrix(c(colLabs,rep("",length(MSEoutNames)-length(colLabs))),ncol=length(colLabs),byrow=TRUE))
  rowLabsAll = as.character(matrix(c(rep("",length(MSEoutNames)-length(rowLabs)),rowLabs),nrow=length(rowLabs),byrow=FALSE))
  
  nrow <- length(rowLabs)
  ncol <- length(colLabs)
  
  par(mfcol=c(nrow,ncol),mgp=mgp,mar=mar,oma=oma,tck=-0.01,xpd=FALSE,xaxs="i",yaxs="i")
  
  ### LOOP THROUGH MSE OBJECTS
  for(outName_i in MSEoutNames){
    topPanels <- (length(scenNames)*(seq_along(OMNames)-1)+1)
    bottomPanels <- length(scenNames)*seq_along(OMNames)
    leftPanels <- seq_along(scenNames)
    rightPanels <- tail(seq_along(MSEoutNames),length(scenNames))
    
    i <- which(MSEoutNames==outName_i)
    out_i <- get(outName_i)
    MPNames <- out_i@MPs
    if(is.null(sims)){
      sims <- 1:out_i@nsim
    }
    if(!is.expression(myvarName)){
      myvar <- slot(out_i,myvarName)  
    }else{
      myvar <- eval(myvarName) 
    }
    proyears_i <- out_i@proyears
    
    legend_panel <- bottomPanels[1]
    
    if(!is.null(myvarDenomExpr)){
      myvarDenom <- eval(myvarDenomExpr)
      myvar <- myvar/myvarDenom # Scale myvar
    }
    
    xaxt_i <- ifelse(i%in%bottomPanels,"n","n")
    yaxt_i <- ifelse(i%in%leftPanels,"n","n")
    
    MPNamesSubix <- match(MPNamesSub,MPNames)
    
    if(plot_type=="huynh"){
      if(!is.null(dim(myvar))){
        myvarsub <- myvar[sims,match(MPNamesSub,MPNames),]
        myvarsubStat <- apply(myvarsub,2,function(x){apply(x,2,FUN)})
      }else{
        myvarsubStat <- matrix(rep(myvar,length(MPNamesSub)),ncol=length(MPNamesSub))  
      }
      dimnames(myvarsubStat) <- list(1:nrow(myvarsubStat),MPNamesSub)
      
      if(!is.null(MPNamesRef)){
        myvarsubStat <- FUNRef(myvarsubStat,MPNamesRef)
      }
      
      if(is.null(ylim_user)){
        ylim <- range(myvarsubStat,na.rm=TRUE)
      }
      if(is.null(xlim_user)){
        xlim <- range(as.numeric(rownames(myvarsubStat)),na.rm=TRUE)
      }
      
      matplot(as.numeric(rownames(myvarsubStat)),myvarsubStat,
              col=cols,lty=ltys,lwd=lwds,type=types,
              xaxt=xaxt_i,yaxt=yaxt_i,
              xlim=xlim,ylim=ylim,xlab="",ylab="")
      #title(main=colLabsAll[i],line=1)
      
      grid(nx=NA,ny=NULL)
      abline(v=assessmentYears,lty="1111",lwd=2)
    }
    
    if(plot_type=="sims"){
      cols <-  ROGBP(length(sims))
      ltys <-  1:length(MPNamesSub)
      pchs <-  NA
      lwds <-  3
      types <- "l"
      
      myvarsub <- myvar[sims,match(MPNamesSub,MPNames),,drop=FALSE]
      dimnames(myvarsub) <- list(sims,MPNamesSub,1:proyears_i)
      
      if(!is.null(MPNamesRef)){
        myvarsub <- FUNRef(myvarsub,MPNamesRef)
      }
      
      if(is.null(ylim_user)){
        ylim <- range(myvarsub,na.rm=TRUE)
      }
      if(is.null(xlim_user)){
        xlim <- c(1,dim(myvarsub)[3])
      }
      
      for(MP_j in MPNamesSub){
      j <- match(MP_j,MPNamesSub)
      myvarsub_j <- myvarsub[,j,]
      lty <- ltys[j]
      if(j==1){
      matplot(t(myvarsub_j), 
              col=cols,lty=lty,lwd=lwds,type=types,
              xaxt=xaxt_i,yaxt=yaxt_i,
              xlim=xlim,ylim=ylim,xlab="",ylab="")
      #title(main=colLabsAll[i],line=1)
      
      grid(nx=NA,ny=NULL)
      abline(v=assessmentYears,lty="1111",lwd=2)
      }else{
        matpoints(t(myvarsub_j), 
                col=cols,lty=lty,lwd=lwds,type=types
                )        
      }
    }
    }
    
    if(plot_type=="addind"){
      if(!is.null(dim(myvar))){
        myvarsub <- myvar
        myvarsubStat <- apply(myvarsub,2,function(x){apply(x,2,FUN)})
      }else{
        # myvarsubStat <- matrix(rep(myvar,length(MPNamesSub)),ncol=length(MPNamesSub))  
      }
      # dimnames(myvarsubStat) <- list(1:nrow(myvarsubStat),MPNamesSub)
      
      # if(!is.null(MPNamesRef)){
      # myvarsubStat <- FUNRef(myvarsubStat,MPNamesRef)
      # }
      
      if(is.null(ylim_user)){
        ylim <- range(myvarsubStat,na.rm=TRUE)
      }
      if(is.null(xlim_user)){
        xlim <- c(1,nrow(myvarsubStat))
      }
      
      matplot(myvarsubStat,
              col=cols,lty=ltys,lwd=lwds,type=types,
              xaxt=xaxt_i,yaxt=yaxt_i,
              xlim=xlim,ylim=ylim,xlab="",ylab="")
      #title(main=colLabsAll[i],line=1)
      
      grid(nx=NA,ny=NULL)
      abline(v=assessmentYears,lty="1111",lwd=2)
    }
    
    if(plot_type=="boot"){
      if(is.null(ylim_user)){
        ylim <- range(myvar,na.rm=TRUE)
      }
      bamExtras::plot_boot_vec(myvar,ylim=ylim)
      abline(v=out_i@nyears,lwd=2,lty=2) # Plot vertical line at end of historical period
    }
    
    if(plot_type=="phase"){
      ltys <-  setNames(rep(0,length(MPNamesSub)),MPNamesSub)
      # pchs <-  setNames(rep(1,length(MPNamesSub)),MPNamesSub)
      lwds <-  setNames(rep(0,length(MPNamesSub)),MPNamesSub)
      types <- setNames(rep("p",length(MPNamesSub)),MPNamesSub)  
      xlabel <- "F/Fmsy"  
      ylabel <- "SSB/SSBmsy"  
      
      if(is.null(proyear_phase)){
        proyear_phase <- out_i@proyears
      }
      Bstatus <- out_i@SB_SBMSY[,MPNamesSubix,proyear_phase]
      Fstatus <- out_i@F_FMSY[,MPNamesSubix,proyear_phase]
      
      matplot(Fstatus,Bstatus,
              col=cols,lty=ltys,lwd=lwds,pch=pchs,type=types,
              xaxt=xaxt_i,yaxt=yaxt_i,
              xlab="",ylab="",xlim=xlim,ylim=ylim
      )
      points(apply(Fstatus,2,median),
             apply(Bstatus,2,median),pch=21,cex=3,lwd=2,col="black",bg=cols)
      
      abline(v=1,lty=3)
      abline(h=1,lty=3)
      
    }
    
    if(plot_type=="tradeoff"){
      out_i <- get(outName_i)
      
      PMx_out <- PMx(out_i)
      PMy_out <- PMy(out_i)
      
      plot(PMx_out@Mean[MPNamesSubix],PMy_out@Mean[MPNamesSubix],pch=16,col=cols,
           xaxt=xaxt_i,yaxt=yaxt_i,
           xlab="",ylab="",
           xlim=xlim,ylim=ylim
      )
      
      # legend("top",legend=MPs,col=MPcols,pch=16,bty="n")
      
      ltys <-  setNames(rep(0,length(MPNamesSub)),MPNamesSub)
      lwds <-  setNames(rep(0,length(MPNamesSub)),MPNamesSub)
      types <- setNames(rep("p",length(MPNamesSub)),MPNamesSub)  
      xlabel <- PMx_out@Caption 
      ylabel <- PMy_out@Caption
      
    }
    
    if(!is.null(hline)){
      abline(h=hline)
    }
    
    if(is.null(x.at_user)){
      x.at <- pretty(par("usr")[1:2])
    }
    if(is.null(y.at_user)){
      y.at <- pretty(par("usr")[3:4])
    }
    
    if(i%in%topPanels){
      col.ct <- which(topPanels==i)
      mtext(colLabsAll[i],side=3,line=1,outer=FALSE,font=2)
    }
    
    if(i%in%rightPanels){
      row.ct <- which(rightPanels==i)
      rowcenter <-par("usr")[3]+diff(par("usr")[3:4])/2
      
      text(x=par("usr")[2]*1.05,y=rowcenter,rowLabsAll[i],
           srt=-90,xpd=NA,font=2)
    }
    
    if(label_axes=="all"){
      x.at.i <- x.at
      axis(side=1,at=x.at.i)
      y.at.i <- y.at
      axis(side=2,at=y.at.i)
    }
    if(label_axes=="leftbottom"){
      if(i%in%bottomPanels){
        x.at.i <- x.at
        if(i>bottomPanels[1]){x.at.i <- x.at[-1]}
        axis(side=1,at=x.at.i,outer=TRUE)
      }
      if(i%in%leftPanels){
        y.at.i <- y.at
        if(i>leftPanels[1]){y.at.i <- y.at[1:(length(y.at)-1)]}
        axis(side=2,at=y.at.i,outer=TRUE)
      }
    }  
    
    if(i==legend_panel&plot_type%in%c("huynh","phase","tradeoff")){
      legend(legend_x,
             inset=legend_inset,
             horiz=TRUE,bty="n",
             legend=MPNamesSub_legtext[legtext_o],
             col=cols[legtext_o],lty=ltys[legtext_o],lwd=lwds[legtext_o],pch=pchs[legtext_o],
             cex=legend_cex,
             xpd=NA,
             x.intersp=0.5
      ) 
      
    }
  }
  
  mtext(ylabel,side=2,line=par("mgp")[1],outer=TRUE)
  mtext(xlabel,side=1,line= par("mgp")[1],outer=TRUE)
  
  