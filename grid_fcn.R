# random number determining location of pixels in each group
randloc_t <- sample.int(nobs, size = nobs, replace = FALSE)
randloc_u <- sample.int(nobs, size = nobs, replace = FALSE)

#subset dataframe based on year and treatment to give treated/untreated pixels for a given year
for (i in 1:(years*2)){
  defor_df$treat <- as.numeric(defor_df$treat)
  #subset to only treated obs in year i
  tyear <- subset(defor_df, variable ==i & treat ==1)
  #attach random location column
  tyear <- cbind(tyear, randloc_t)
  #sort based on random location column
  tyear <- tyear[order(randloc_t),]
  
  uyear <- subset(defor_df, variable ==i & treat ==0)
  uyear <- cbind(uyear, randloc_u)
  uyear <- uyear[order(randloc_u),]
  
  #splits into "counties" based on desired size of county x
  treatgrid <- split(tyear, ceiling(seq_along(tyear$idx)/x))
  untreatgrid <- split(uyear, ceiling(seq_along(uyear$idx)/x))
  
  #  for j in 1:length(treatgrid){
  #   deforrate_t <- mean(treatgrid[j]) 
  #  deforrate_u <- mean(untreatgrid[j])
  #  }
  
  
  #name based on year
  t_name <- paste("treatgrid", i, sep = "_")
  u_name <- paste("untreatgrid", i, sep = "_")
  assign(t_name, treatgrid)
  assign(u_name, untreatgrid)
  rm(tyear, uyear, treatgrid, untreatgrid)
  
}