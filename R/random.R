
library(RColorBrewer)
rf <- colorRampPalette(rev(brewer.pal(11,'Spectral')))
r <- rf(32)
library(hexbin)
h <- hexbin(df <- data.frame(foo[i,3], -Mk[i]))
plot(h)
plot(h, colramp=rf)

