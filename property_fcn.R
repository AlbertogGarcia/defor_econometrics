# function to aggregate pixels to property level, adding in random perterbations at property level



# proposed layout of function:
# 1) create polygon with size of landscape
# 2) generate random points in polygon and perform voneroi tesselation to determine county boundaries
# 3) randomly assign counties to treatment or control status
# 4) randomly assign units to location on landscape depending on treatment status
# 5) calculate average defor rate within each county
# 6) 

# #now within the landscape, we generate 10 points from which to generate our voronoi tesselations
# points = 10
# 
# x <- sample(1:rootn,points)
# y <- sample(1:rootn,points)
# 
# pol <- data.frame(x, y)
# ggplot(pol,aes(x,y)) +
#   stat_voronoi(geom="path") +
#   geom_point()















# ending function
# }



# # random number determining location of pixels in each group
# randloc_t <- sample.int(nobs, size = nobs, replace = FALSE)
# randloc_u <- sample.int(nobs, size = nobs, replace = FALSE)
# 
# #subset dataframe based on year and treatment to give treated/untreated pixels for a given year
# for (i in 1:(years*2)){
#   defor_df$treat <- as.numeric(defor_df$treat)
#   #subset to only treated obs in year i
#   tyear <- subset(defor_df, year ==i & treat ==1)
#   #attach random location column
#   tyear <- cbind(tyear, randloc_t)
#   #sort based on random location column
#   tyear <- tyear[order(randloc_t),]
#   
#   uyear <- subset(defor_df, year ==i & treat ==0)
#   uyear <- cbind(uyear, randloc_u)
#   uyear <- uyear[order(randloc_u),]
#   
#   x = 1000*d/30
#   
#  
#   #name based on year
#   #t_name <- paste("", i, sep = "_")
#   #u_name <- paste("", i, sep = "_")
#   assign(t_name, treatgrid)
#   assign(u_name, untreatgrid)
#   rm(tyear, uyear, treatgrid, untreatgrid)
#   
# }
# 

