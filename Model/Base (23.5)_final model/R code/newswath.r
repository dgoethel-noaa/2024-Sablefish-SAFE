newswath<-function(sims=mcmc/1000,styear,endyear,ylab="Spawning biomass (kt)",B40plot=1)
{
    #Function to plot quantiles of time-series  
    #Arguments to function are:
    # 1) object
    # 2) number of sim trajectories to overlay
    # 3) yaxis dim                               ,48,30,yr.first=1960
    # 4) label for y axis                        
    # 5) relative to first year                  
    # 6) first year                  
    # 7) last year   
    # 8) B40plot, setting for plotting lines of B40 and B35 by species
    #    0 = no lines for B40 or B35 on swathplat, 
    #    1 = plot lines based on species in Data2R.r fnc               
    #sims2<-sims
    # sims2<-matrix(as.numeric(unlist(sims[(9+2*ys-2+nages):(9+3*ys-1+13+nages)])),4501,ys+15)
    nyrs<-endyear-styear+1
    recyears<-nyrs+28 # 28 prerecruit variables
    yr.first<-styear
    yr.last<-endyear+15
    sims2<-matrix(as.numeric(unlist(sims$mcmc[1:length(sims$mcmc)])),4501,length(sims$mcmc))
    sims2<-cbind(sims2[,(nyrs+14+1):(2*nyrs+14)],sims2[,(3*nyrs+14+1):(3*nyrs+2*14+1)]) # both spbio and bpbio_proj move forward by +1
    .Options$warn <- -1  # suppress annoying warning messages
    
    ssb.quantiles<- apply(sims2,2,function(x){quantile(x,seq(.025,.975,.05))})
    #ssb.quantiles<- apply(sims2,2,function(x){quantile(x,seq(.1,.90,.10))})
    print(summary(ssb.quantiles))
    n <- nrow(sims2)
    nquant=nrow(ssb.quantiles)
    ssb.med2004 <- 1.
    
    #win.graph(height=8,width=10)
    #par(mai=c(0,1,0,.5),omi=c(1.2,0,1,0))
    par(cex.axis=1.2,cex.lab=1.3)
    
    #  pall<-c( "lightblue2","lightblue1","lightblue",  "skyblue1" ,  "skyblue2", "skyblue3" , 
    #   "deepskyblue1","deepskyblue2" , "darkblue","deepskyblue3" )
    #    col=paste("gray",round(seq(90,5,-90*2/nquant))[i],sep=""),border=NA)
    
    #Options for color ramp using vcd library,   
    #pall <- rainbow(9)
    pall <- rainbow(9,s=0.65,v=0.85,start=0.7,end=0.95)
    #pall <- rainbow(9,start = 0.7, end = 0.83)
    #pall <- heat.colors(9)
    #pall <- terrain.colors(9)
    #pall <- topo.colors (9)
    #pall <- cm.colors (9)	  
   

    ylim <- c(0,max(ssb.quantiles/ssb.med2004))
    plot(yr.first:yr.last,ssb.quantiles[(nquant+1)/2,]/ssb.med2004,type="l",ylim=ylim,
         xlab="Year",ylab=ylab, xlim=c(yr.first-.5,yr.last+.5),cex=1,yaxs="i")
    for(i in 1:((nquant-1)/2)){
        polygon(c(yr.first:yr.last,yr.last:yr.first),c(ssb.quantiles[i,],
                                                       rev(ssb.quantiles[nquant-i+1,]))/ssb.med2004,
                col=paste(pall[i],sep=""),border=NA)
        
    } 
    lines(yr.first:yr.last,ssb.quantiles[(nquant+1)/2,]/ssb.med2004,type="l",lwd=2,col="white")
    lines(yr.first:yr.last,colMeans(sims2)/ssb.med2004,type="l",lwd=2,col="green")
    
    if(B40plot>0) {
        lines(yr.first:yr.last,	rep(B40,length(yr.first:yr.last)), type="l", lty=2, lwd=3, col="red")	
        lines(yr.first:yr.last, rep(B35,length(yr.first:yr.last)), type="l", lty=1, lwd=3)
        text(1965,B40+32,expression(paste(italic(B)["40%"])),col="red",cex=1.6)
        text(1965,B35-30,expression(paste(italic(B)["35%"])),col="black",cex=1.6)
        text(yr.last-15,B40+24,as.character(endyear+1),col="dark green",cex=1.3)
    }
}
