
# Function to plot color bar
color.bar <- function(lut, min, max=-min, nticks=11, ticks=seq(min, max, len=nticks), title='') {
  scale = (length(lut)-1)/(max-min)

  # dev.new(width=1.75, height=5)
  plot(c(0,10), c(min,max), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='')
  title(title, line = 1,cex.main = 1)
  axis(2, ticks, las=1,line = -4)
  for (i in 1:(length(lut)-1)) {
    y = (i-1)/scale + min
    rect(6,y,8,y+1/scale, col=lut[i], border=NA)
  }
}
