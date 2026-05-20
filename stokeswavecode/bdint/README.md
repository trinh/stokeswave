qt_meshsurface_phil
Sept 2015
Modification of qt_mesh_surface to do both top and lower surfaces

### 18 May 2026

Extra notes added due to Henry asking about the correctedness of the series vs. boundary integral solver. 

* Looking at these old codes, I see two lines. The `bdint_` codes are likely written at a previous state, and these seem to solve the analytic continuation in z(f). At some point, this was likely converted to the (q, theta) formulation, and so this is the `qt_bdint_` set of codes

* I can see that occasioanlly some of the `boundaryintegral` scripts are broken. For example, `boundaryintegral_single` references a `delta(:)` vector 

* In many scripts, the data file location should be changed to e.g. `F = load('./../../stokeswavedata/series_N100_ep0.2.mat');`

* Would first suggest running `qt_boundarycomparison.m` which plots the series data and then compares with a boundary integral solution along the real axis, showing excellent fit. 

* Then run `qt_bdint_path_call_henry_20260518` for a check of the complex boundary integral solver




