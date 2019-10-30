# PCPR2 on either all coffee samples (n=174) or all coffee samples except instants (n=150)
# Final analysis for coffee brew paper

# Example from pcpr2 package
library(devtools)
install_github("JoeRothwell/pcpr2")

library(pcpr2)
output <- runPCPR2(transcripts, Y_metadata)
plot(output)

# Run PCPR2_coffee.R lines 7-38
pcpr2.coffee <- runPCPR2(X_DataMatrixScaled, Z_InterestFactors)
plot(pcpr2.coffee)
summary(pcpr2.coffee)

# plotProp is now deprecated
plotProp(pcpr2.coffee)
