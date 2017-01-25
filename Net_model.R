Net.model <- function(e=edges,farms=attr1,mu=.5,beta=.1,alph=2){  
   #mu=Probability that each animal in the batch is infected
   #alph = probability of local spread, determines decay over distance
   e$Date <- as.Date(e$Date)
   names(e)[1:2]<-c("origin","dest")
   
   #get list of farms for the network (including isolates)
   farms$name <- farms$farm.id
   
   t <- seq.Date(from=min(e$Date),to=max(e$Date),by=1)
   timestep <- length(t)
   
   #file to keep track of each farms state.  Updates state each timestep.  Starts with 1 randomly infected farm
   farms$state <- sample(c(rep(0,nrow(farms)-1),1),nrow(farms))#0: Uninfect, 1: infected
   farms$time.I <- ifelse(farms$state==1,0,NA)  #for recording time infected.  NA for all except initially infected farm
   farms$tmode.local <- "NA"  #Record if local transmission
   farms$tmode.move <- "NA"   #Record if movement transmission
   farms$mode <- "NA"       #Record if movement & local transmission
   
   #file for keeping track of S and I numbers
   track <- data.frame(time=min(t),S=c(nrow(farms)-1),I=c(1))
   S <- nrow(farms)-1   #object to keep track of how numbers are changing each timestep
   I <- 1                #object to keep track of how numbers are changing each timestep
   
   
   rbinom.v <- function(x){rbinom(1,1,x)}
   
   plot(farms$long[farms$state==0],farms$lat[farms$state==0],xlab="Longitude",ylab="Lattitude")
   points(farms$long[farms$time.I==0],farms$lat[farms$time.I==0],col="slateblue",pch=19,cex=1.3)
   legend("topright",legend=c("Uninfected","Infected through movements","Infected through local spread","Index case"),col=c("black","red","green4","slateblue"),pch=c(1,17,19,19),cex=.8,bty="n")
   
   
   for(i in 2:timestep){
      #create new row for next timestep with same values as previous (will be updated later if transmission occurs)
      row <- data.frame(time=t[i],S=S,I=I)
      track <- rbind(track,row)
      
      points(farms$long[farms$tmode.move=="Movement"],farms$lat[farms$tmode.move=="Movement"],col="red",pch=17,cex=1.3)
      points(farms$long[farms$tmode.local=="Local"],farms$lat[farms$tmode.local=="Local"],col="green4",pch=19,cex=1.3)
      
      ###LOCAL SPATIAL SPREAD
      if(alph>0) {#only run if alph set to  >0
         f.i <- subset(farms,state==1,select=c("farm.id","lat","long")) #get list of infected farms
         f.u <- farms[farms$state==0,1]  #get list of uninfected farms
         f.m <- merge(f.i,f.u)  #merge infected farms with all uninfected farms
         #add Uninfected farms UTM
         f.m$long.u <- farms$long[match(f.m$y,farms$farm.id)]
         f.m$lat.u <- farms$lat[match(f.m$y,farms$farm.id)]
         f.m$dist.sq <- (f.m$long-f.m$long.u)^2 + (f.m$lat-f.m$lat.u)^2 #measures as d^2 because quicker
         f.r <- subset(f.m,dist.sq<20000^2) # filter for only farms that are within 20 km
         f.r$km <- ((f.r$dist.sq)^.5)/1000 # get distance in km
         f.r$y <- as.character(f.r$y)
         #It is farm y that may become infected by all the farms in name col of f.r
         f.r$prob <- 1-mu*beta*exp(-alph*(f.r$km+.01)) #Prob that trans does NOT occur from farm.x to y, add +.01 so that any farms with distance 0 are not 100% likely to transmit
         #multiply by mu to rescale the transmission risk based on within-farm prev.
         if(nrow(f.r)>0){ #skip if no farms within specified distance
            trans.s <- data.frame(PROB.cum =tapply(f.r$prob,f.r$y,prod))#multiply all prob together for each uninfected farm
            trans.s$name.y <- row.names(trans.s)
            trans.s$m <- 1-trans.s$PROB.cum
            trans.s$trans <- c()
         }}
      #will compute transission later because we don't want farms to become infectious until next timesetp (we still need to compute trans from movement)
      
      
      ###MOVEMENTS
      M.o <- nrow(e[e$Date==t[i],]) #observed number of movements on time t
      
      if(M.o>0){  #skips to below if no movements occur
         edges.t <- subset(e,Date==t[i])
         
         #bring in state of senders and receivers
         edges.t$state1 <- farms$state[match(edges.t$origin,farms$name)]
         edges.t$state2 <- farms$state[match(edges.t$dest,farms$name)]
         if(sum(edges.t$state1,na.rm=T)>0){  #skip to below if no infected farms moved any animals
            #filter edgelist so only includes movements from I --> S farms
            edges.f <- subset(edges.t,state1==1 & state2==0)
            
            if(nrow(edges.f)>0){ #skip to below if no I-->S movements
               #calculate state transitions
               edges.f$PROB <- (1-mu)^edges.f$batch.size
               trans <- data.frame(PROB.cum =tapply(edges.f$PROB,edges.f$dest,prod))#multiply all prob together
               trans$name <- row.names(trans)
               trans$m <- 1-trans$PROB.cum
               #get transition
               trans$trans <- sapply(trans$m,rbinom.v)
               
               #if transmition occurs, update state in farms
               trans.names <- trans$name[trans$trans==1]
               farms$state <- ifelse(farms$name %in% trans.names,1,farms$state)
               farms$tmode.move <- ifelse(farms$name %in% trans.names,"Movement",farms$tmode.move)
               farms$time.I <- ifelse(farms$name %in% trans.names,i,farms$time.I)
               
               if(sum(trans$trans)>0){#update farm counts
                  S <- S - sum(trans$trans) #decrease S count by one farm
                  I <- I + sum(trans$trans) #increase I count by one farm
               } 
            }}}
      
      #now update farms from the local spatial transmisison
      if(alph>0) {#only run if alph>0
         if(nrow(f.r)>0){ #skip if no farms within specified distance
            #get transition
            trans.s$trans <- sapply(trans.s$m,rbinom.v)
            
            #if transmition occurs, update state in farms
            trans.s.names <- trans.s$name[trans.s$trans==1]
            farms$state <- ifelse(farms$name %in% trans.s.names,1,farms$state)
            farms$tmode.local <- ifelse(farms$name %in% trans.s.names,"Local",farms$tmode.local)
            farms$time.I <- ifelse(farms$name %in% trans.s.names,i,farms$time.I)
            
            if(sum(trans.s$trans)>0){#if transmition occurs, update state in farms
               S <- S - sum(trans.s$trans) #decrease S count by one farm
               I <- I + sum(trans.s$trans) #increase I count by one farm
            }
         }}
      #Update counts of tracker   
      track$S[i] <- S
      track$I[i] <- I   
   }
   
   #Record mode of transmission for farms in a single column.
   farms$mode <- ifelse(farms$tmode.local=="Local" & is.na(farms$tmode.move)==F,"Local",ifelse(is.na(farms$tmode.local)==F & farms$tmode.move=="Movement","Movement",ifelse(farms$tmode.local=="Local" & farms$tmode.move=="Movement","Movement and Local",NA)))
   
   return(list(farms,track))
}