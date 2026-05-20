2026-05-20: The dual mesh for both surface 0 and surface A look OK. 

The trick to finding the diagonal singularities is as follows. 

Look into stokeswavepub/contourep0p5 to see the meshes that were used in the publication Fig. 10

You will see that diagonal singularities appears in the file 'dualmesh.AA_.fig' These sigularities were around phi = 0.5, psi = -0.4. 

The problem is that in finite depth, the mirrored boundaries have depth 

Q = -log(r0)/(2*pi) or r0 = exp(-2 pi Q)

So if a singularity lies at Im(f) = -0.4, this corresponds to a value of r0 of around 0.08. In other words, unless you extend the code to allow you to see beyond the thresholds -Q < Im(f) < Q, you need to take r0 very small to see the diagonal singularities (unless there is more interesting stuff that we have not anticipated). 