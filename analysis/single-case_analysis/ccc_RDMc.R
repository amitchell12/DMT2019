#define function for case control comparisons
#xmean = control x mean, xsd = control x sd, ymean = control y mean, ysd = control y sd, n = control sample size, x = case x score, y = case y score_
ccc <- function(xmean, xsd, ymean, ysd, rxy, n, x, y){
  
  #z scores of case
  zx <- (x - xmean)/xsd
  zy <- (y - ymean)/ysd
  
  #one-tailed tests of deficit on tasks x and y
  TD_x <- zx/sqrt((n+1)/n) <= qt(p = .05, df = n-1)
  TD_y <- zy/sqrt((n+1)/n) <= qt(p = .05, df = n-1)
  
  #TESTS OF DIFFERENCE USING UDT (note commented out alternative to use RSDT)
  
  #UDT: Crawford, J. R., & Garthwaite, P. H.  (2005).
  #Testing for suspected impairments and dissociations in single-case studies in neuropsychology:
  #Evaluation of alternatives using Monte Carlo simulations and revised tests for dissociations.
  #Neuropsychology, 19, 318-331.
  #UDT <- ((x - xmean) - (y - ymean))/sqrt(((xsd*xsd) + (ysd*ysd) - (2*xsd*ysd*rxy))*((n+1)/n))
  #UDT_xy <- abs(UDT) >= qt(p = .025, df = n-1, lower.tail = FALSE) #two-tailed UDT
  #UDT_2p <- 2*pt(abs(UDT), df = n-1, lower.tail = FALSE)
  
  #RSDT: Crawford, J. R., & Garthwaite, P. H.  (2005).
  #Testing for suspected impairments and dissociations in single-case studies in neuropsychology:
  #Evaluation of alternatives using Monte Carlo simulations and revised tests for dissociations.
  #Neuropsychology, 19, 318-331.
  a <- (1 + rxy)*(1 - (rxy*rxy))
  b <- (1-rxy)*(4*((n-1)*(n-1)) + 4*((1+rxy)*(n-1)) + (1+rxy)*(5+rxy))
  c <- (-2*((zx - zy)*(zx - zy))) * ((n*((n-1)*(n-1)))/(n+1))
  y <- ((-b + sqrt((b*b) - (4*a*c)))/(2*a))^.5
  RSDT_xy <-  y >= qt(p = .025, df = n-1, lower.tail = FALSE) #two-tailed RSDT lower
  RSDT_2p <- 2*pt(y, df = n-1, lower.tail = FALSE)
  
  #Z_DCC (effect size of Difference Case-Controls): Crawford, J. R., Garthwaite, P. H., and Porter, S. (2010).
  #Point and interval estimates of effect sizes for the case-controls design in neuropsychology:
  #Rationale, methods, implementations, and proposed reporting standards.
  #Cognitive Neuropsychology, 27, 245-260. 
  #note, negative Z_DCC indicate relative impairment on X, positive indicate relative impairment on Y
  Z_Dcc <- (zx - zy)/sqrt(2-(2*rxy))
  
  #specify output vector
  ccc_out <- list(Zcc_x = zx, Zcc_y = zy, Z_Dcc = Z_Dcc, TD_x = TD_x, TD_y = TD_y, RSDT_xy = RSDT_xy,  RSDT_2p = RSDT_2p)
  return(ccc_out)
}