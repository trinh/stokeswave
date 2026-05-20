## Performs the analytic continuation using the boundary integral

**Original writing:** sometime in 2017

### 2026 re-write

We rewrote scripts in this directory. The main changes were as follows.


 - We now have `qt_finite_boundarycomparison.m` does computes both series and boundary-integral solutions on the real axis and verifies they agree
 - The `qt_bdint_path_call.m` script attempts to solve for the complexification of Omega along a desired path in the complex plane. 

Labeled everything with _finite please to avoid confusion!

The other thing that is done is to move all functions into a central functions folder. This allows the creation of subfolders for working individual studies, rather than over-writing scripts. 