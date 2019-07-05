library(DeclareDesign)
  
DIDbias <- vector()

for (p in 1:100){
  
  
  b0 <- .3
  b1 <- .2
  b2 <- -.1
  b3 <- -.23
  
  sd_a <- 0.5
  sd_v <- 0.5
  
  ATT <- pnorm(b0+b1+b2+b3,0, 1 ) - pnorm(b0+b1+b2,0, 1 )
  
  panels <- fabricate(
    pixels = add_level(N = nobs, a_i = rnorm(N, 0, sd_a), treat = rbinom(N, 1, 0.5)),
    year = add_level(N = (years*2), nest = FALSE),
    obs = cross_levels(
      by = join(pixels, year),
      post = ifelse(year > years, 1, 0),
      v_it = rnorm(N, 0, std_v),
      ystar = b0 + b1*treat + b2*post + b3*treat*post + a_i + v_it,
      y_it = ifelse(ystar > 0, 1, 0)
    )
  )
  
  
  
   probitcoeff  <- glm(y_it ~  post*treat,
                       family = binomial(link = "probit"),
                       data = panels
   )$coefficients

   tr <- sum(probitcoeff)
   untr <- probitcoeff[1]+probitcoeff[2]+probitcoeff[3]
   probt <- pnorm(tr)
   probu <- pnorm(untr)
   coeff1[p] <- probt-probu
  
  DIDbias[p] <- lm(y_it ~ post*treat, data = panels)$coefficients[4]-ATT
  
}

summary(DIDbias)


