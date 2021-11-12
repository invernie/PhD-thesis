# Checking Figure 3 in Franks and Deneubourg (1997)

set.seed(0)

s0 <- 50  # n of stones being carried at time 0
p <- 0.05 # prob of encountering a stone
d <- 0.046 # prob of deposition

niter <- 100
rsqv <- rep(0,niter)  # vector of r squared of the linear models

for (iter in 1:niter) {

  antv <- rep(1,s0) # vector recording what ants are still carrying
  stonev <- rep(0,s0) # vector recording number of encountered stones
  s <- s0
  
  
  while (s >0) {  # until all stones have been deposited
    
    for (i in 1:s0){  # go through each ant that was carrying to start with
      
      if (antv [i] > 0) {  # if it is still carrying
        
        if (runif(1)<= d) {  # check if it deposits
          
          s <- s-1
          antv [i] <- 0
          
        } else {
          
          if (runif(1) <= p) {  # check if it encounters a stone
            
            stonev [i] <- stonev [i] + 1 
            
          }
          
        }
        
      }
      
    }
    
  }
  
  
  br <- seq(0,max(stonev),1)
  h <- hist(stonev, breaks = br)
  sdep <- h$counts
  scarr <- c()
  for (x in 1:length(sdep)) {
    z <- 0
    for (xx in 1:x) {
      z <- z+sdep[xx]  # count how many stones have been deposited at previous values of n stones encountered
    }
    scarr <- c(scarr, s0-z) # subtract from starting n of stones to calculate how many are remaining
  }
    
  m <- lm(c(log(scarr[2:(length(scarr)-1)]),0)~br[2:(length(sdep))]) # exclude zero stone data point(as in F&D)
      
  rsqv[iter] <- summary(m)$r.squared
  
}

# plot last sim
dev.new(width=5, height=5)
par(mar = c(6,6,3,3), mfrow = c(1,1))
windowsFonts(A = windowsFont("Times New Roman"))
#plot(y = scarr,  x = br[1:(length(sdep))], ylab = "N of stones still carried", xlab = "N of stones encountered")
plot(y = log(scarr),  x = br[1:(length(sdep))], xlim = c(-0.2, 4.5), ylim = c(-0.2,3), family = 'A', ylab = " ", xlab = " ", cex.lab = 1.5, cex.axis = 1.5, tck = 0.04, pch = 20, cex = 2)
clip(-0.1,4.2,-0.05,3)
abline(m, lwd = 2)

mean(rsqv)
median(rsqv)
sd(rsqv)

length(which(rsqv >= 0.92))

