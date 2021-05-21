full_summary <- read.csv("full_summary.csv")#, check.names = FALSE)
# Looks better when there is an outer margins

par(oma=c(1,0,1,1))

labels <- list(#"Model:" = " ",
  "Unit of analysis:" = c("pixel", "grid", "property", "county"),
  "Fixed effects:" = c("grid FE", "property FE", "county FE", "treatment FE"),
  "SE structure:" = c("clustered at pixel", "clustered at grid", "clustered at property", "clustered at county"))

full_summary <- subset(full_summary, pixel.fe!="TRUE" |is.na(pixel.fe))


full_summary <- rbind(full_summary[1:4,], full_summary[8:10,], full_summary[5:7,], 
                      full_summary[11:15,], full_summary[19:21,], full_summary[16:18,],
                      full_summary[22:26,], full_summary[30:32,], full_summary[27:29,],
                      full_summary[33:37,], full_summary[41:43,], full_summary[38:40,],
                      full_summary[44,])



f <- function(x) gsub("^(\\s*[+|-]?)0\\.", "\\1.", as.character(x))

coverage <- full_summary$cover
coverage[is.na(coverage)] <- 0
c_print <- f(round(full_summary$cover, digits = 2))
c_print[is.na(c_print)] <- " "

RMSE <- full_summary$RMSE
RMSE[is.na(RMSE)] <- 0
RMSE <- as.numeric(RMSE)

RMSE_print<- f(round(full_summary$RMSE, digits =3))
RMSE_print[is.na(RMSE_print)] <- " "

column_indic <- c(1:10, " ", 11:20, " ", 21:30, " ", 31:40, " ")

test <- subset(full_summary, select=-c(X,pixel.fe, cover, RMSE, sigma_p))%>%
  filter(pixel=="TRUE"|pixel=="FALSE")

index.ci <- match(c("q05","q95"), names(test))

# One could also add information about model fit to this chart

topline = -0.035
midline = topline-0.03
ylim <- c(midline-0.015,0.065)
#Create the plot

schart(test,labels, ylim = ylim, index.ci=index.ci, ylab="Bias"#col.est=c("grey50", "blue") #, highlight=c(1:7, 15:21)
       #,band.ref=c(-.05, .04)
       , axes = FALSE, n=10
       #, col.band.ref="#c7e9f9"
) # make some room at the bottom
Axis(side=2, at = 0, labels=TRUE)
abline(h=topline)
abline(h=midline)
abline(v=11, lty="dashed")
abline(v=22, lty="dashed")
abline(v=33, lty="dashed")
lapply(1:length(RMSE), function(i) {
  rect(xleft=i-.3, ybottom=midline, xright=i+.3, ytop=midline+RMSE[i]*1.5, border=NA, col="#E69F00")
  #mtext(paste0(column_indic[i]), side=1, at = i, font=2, cex=.9)#, line=1, at=-1)
  text(x= i, y=midline+RMSE[i]*1.5+0.004, paste0(RMSE_print[i]), col="black", font=1, cex=.7)
  text(x= i, y=midline-0.01, paste0(c_print[i]), col="black", font=1, cex=.825 )
})
# text(x=5#mean(1:nrow(test))
#      , y=topline-0.0075, "RMSE", col="black", font=2)
mtext("RMSE", side=2, at = midline+0.015, font=2, las=1, line=.5)
# text(x=5#mean(1:nrow(test))
#      , y=midline-0.0075, "coverage probability", col="black", font=2)
mtext("Coverage\nprobability", side=2, at = midline-0.01, font=2, las=1, line=.5)
text(x=1
     , y=.065, expression(paste(sigma[p],"=0.0")), col="black", font=2)
text(x=12.5
     , y=.065, expression(paste(sigma[p],"=0.1")), col="black", font=2)
text(x=23.5
     , y=.065, expression(paste(sigma[p],"=0.2")), col="black", font=2)
text(x=34.5
     , y=.065, expression(paste(sigma[p],"=0.3")), col="black", font=2)
legend(x=29.5, y=0.065, col = c("grey60"), legend = c("0.05 to 0.95 quantile\nof estimate distribution"), inset = 0.005,  box.lty=0, cex=0.75,  seg.len=0.37, lty = 1, horiz=TRUE, lwd = 4, bg="transparent")
