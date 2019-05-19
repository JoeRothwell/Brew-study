# PCPR2 on either all coffee samples (n=174) or all coffee samples except instants (n=150)
# Final analysis for coffee brew paper

# Example from pcpr2 package
library(devtools)
install_github("JoeRothwell/pcpr2")

library(pcpr2)
output <- runPCPR2(transcripts, Y_metadata)
plotProp(output)

# Run PCPR2_coffee.R lines 11-38
output.cof <- runPCPR2(X_DataMatrixScaled, Z_InterestFactors)
output.cof

plotProp(output.cof)
