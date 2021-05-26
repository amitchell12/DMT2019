## A.G.Mitchell 26.05.2021
## SUPPLEMENTARY RAD ANALYSIS - X,Y coordinates at peak speed
## is there a difference between groups?

# subtracting cal from reach endpoint
res$LANDx <- res$px - res$calx
res$LANDy <- res$py - res$caly


# plotting to get a look at data
res$POSITION <- factor(res$POSITION)
ggplot(res) + geom_point(aes(x = calx, y = caly, colour = POSITION), shape = 3) +
  geom_point(aes(x = mx, y = my, colour = POSITION)) +
  facet_wrap(~PPT*VIEW)
ggsave('radial-reach_Err.png', plot = last_plot(), device = NULL, 
       path = anaPath, scale = 1, width = 15, height = 10, units = 'in')

## data transformations and calculations
# tranforming to degrees
res$LANDx_deg <- visAngle(size= res$LANDx, distance = 500) #using visual angle function above
res$LANDy_deg <- visAngle(size= res$LANDy, distance = 500)

# absolute error - in mm
res$AE <- sqrt(res$LANDx^2 + res$LANDy^2) #mm

## angular error and amplitude error
# calculating target and end-point position relative to start-point
res$rX <- res$mx - res$sX
res$rY <- res$my - res$sY
res$tX <- res$calx - res$sX
res$tY <- res$caly - res$sY

#recode response as target-relative ERRORS in polar coordinates
res$tANG <- (atan(res$tX/res$tY))*(180/pi)
res$tAMP <- sqrt(res$tX^2 + res$tY^2)
res$rANG <- (atan(res$rX/res$rY))*(180/pi)
res$rAMP <- sqrt(res$rX^2 + res$rY^2)

## calculating angular error and amplitude error
res$ANG_ERR <- res$rANG - res$tANG
res$AMP_ERR <- res$rAMP - res$tAMP

# removing
res <- res[res$EYE_MOVE == 0, c(1:7,17:18,8:15,22:42)]

#save compiled data-set
write.csv(res, "radial-reaching_compiled.csv", row.names = FALSE)
