summary_df <- read.csv("summary_df.csv")#, check.names = FALSE)
# Looks better when there is an outer margins

summary_df <- rbind(summary_df[1:5,], summary_df[9:11,], summary_df[6:8,])

par(oma=c(1,0,1,1))

labels <- list("Unit of analysis:" = c("pixel", "grid", "property", "county"),
               "Fixed effects:" = c("pixel FE", "grid FE", "property FE", "county FE", "treatment FE"),
               "SE structure:" = c("clustered at pixel", "clustered at grid", "clustered at property", "clustered at county"))
f <- function(x) gsub("^(\\s*[+|-]?)0\\.", "\\1.", as.character(x))
coverage <- f(round(summary_df$cover, digits=3))
RMSE <- f(round(summary_df$RMSE, digits=5))
test <- subset(summary_df, select=-c(X, cover, RMSE))

index.ci <- match(c("q05","q95"), names(test))

# One could also add information about model fit to this chart
ylim <- c(-0.035,0.04)

topline = -0.015

midline = topline-0.02-.0025

ylim <- c(midline-.02,0.04)
#bottomline = (min(ylim)+topline)/2
#Create the plot


schart(test,labels, ylim = ylim, index.ci=index.ci, col.est=c("grey50","#D55E00"), ylab="Bias", highlight=1
       ,band.ref=c(-.05, .04)
       , axes = FALSE
       #, col.band.ref="#c7e9f9"
) # make some room at the bottom
Axis(side=2, at = c( 0, 0.03), labels=TRUE)
abline(h=topline)
abline(h=midline)
#abline(h=bottomline)
lapply(1:length(RMSE), function(i) {
 # text(x= i, y=min(ylim)+.002, paste0(i), col="black", font=2, cex = 1)
  mtext(paste0(i), side=1, at = i, font=2, cex=.95)#, line=1, at=-1)
  text(x= i, y=midline+0.005, paste0(RMSE[i]), col="black", font=1, cex=.9)
  text(x= i, y=min(ylim)+0.0025, paste0(coverage[i]), col="black", font=1, cex=.9 )
})
# mtext('coverage\nprobability', side=2, las=1, at = loc_1, adj = 1, font=2, line=1)
# mtext('RMSE', side=2, las=1, at = loc_2, font=2, line=1)
text(x=mean(1:nrow(test))
     , y=midline-.005, "coverage probability (95% desired)", col="black", font=2)
text(x=mean(1:nrow(test))
    , y=topline-.005, "RMSE", col="black", font=2)
legend(x=6.5, y=0.05, col = "grey50", legend = "0.05 to 0.95 quantile \n of estimate distribution", seg.len=0.65, inset = 0.005,  box.lty=0, cex=0.8, lty = 1, lwd = 4, bg="transparent")
