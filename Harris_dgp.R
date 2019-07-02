#### Delgado-esque dgp

nobs <- 10000 #number of obs in each year
years <- 2 #number of years in each period

b0 <- -.4 

b1 <- .75 # trend between periods

b2 <- .3 #pre treatment difference between groups  

ATT <- -.5

b4 <- .25

df <- data.frame(rbinom(nobs, 1, .5))  #assigning treatment with probability .5
names(df)[1] <- paste("G_i")
df <- tibble::rownames_to_column(df)
df <- unite(df, idx, G_i, rowname, sep = "_", remove = FALSE)
df <- subset(df, select = -c(rowname))

df <- do.call("rbind", replicate(years*2, df, simplify = FALSE))

trend <- rep(rnorm(1, 0, 1), nobs)
a <- 1.02


year <- vector()
for(i in 1:(years*2)){
  
  trend <- c(trend, rep(a*trend[i*nobs]+rnorm(1,0,(a-1)^2), nobs))
  year <- c(year, rep(i, nobs))
  
}
trend <- trend[1:(nobs*years*2)]
df$year <- year
df$trend <- trend
df$post <- (df$year > years)*1

std_a <- .1
std_v <- .25

df$a_i <- rep(rnorm(nobs, 0, std_a), years*2) #individual specific
df$v_it <- rnorm(nobs*years*2, 0, std_v)  


df$ystar <- b0 + b1 * df$trend +b2 * df$G_i + ATT * df$G_i* df$post + df$a_i +df$v_it
df$ystar <- b0 + b1 * df$post +b2 * df$G_i + ATT * df$G_i* df$post + df$a_i +df$v_it + b4 *trend
#df$ystar <- b0 + b1 * df$post +b2 * df$G_i + ATT * df$G_i* df$post + df$a_i +df$v_it


df$y_it <- (df$ystar > 0)*1
defor_df <- df


ones <- subset(defor_df, y_it == 1)
ones <-ones[with(ones, order(year)), ]
first <- ones[match(unique(ones$idx), ones$idx),]
first <- subset(first, select = c(idx, year))
names(first)[2] <- paste("defor_year")

defor_df <- merge(defor_df, first, by = "idx", all = TRUE)
defor_df$defor_year[is.na(defor_df$defor_year)] <- (years*2+1)
defor_df$indic <- (defor_df$year - defor_df$defor_year)
defor_df$y_it <- ifelse(defor_df$indic > 0 , NA, defor_df$y_it)
defor_df <- subset(defor_df, select = -c(indic, a_i, v_it, ystar))


defor_df <- defor_df[order(defor_df$idx, defor_df$year),]


lm(y_it ~  post*G_i, 
   data = defor_df
)$coefficients
