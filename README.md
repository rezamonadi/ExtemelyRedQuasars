## Precise Selection of Extremely Red Quasars

### Our final sample of Type1 Boxy CIV lines Extremely Red quasars (T1BERQs) is available [here](https://github.com/rezamonadi/ExtremelyRedQuasars/raw/main/T1BERQSample.fits) as a fits table

### Sample selection 
Sampling from the [emission line catalog](https://datadryad.org/stash/dataset/doi:10.6086/D1H59V) in H17: 
```erqSampling.ipynb```

### Luminosity Matching
Matching the luminosity of the selected sample with the Luminosity of ERQs:
```lum-match.ipynb```
### Sample Statistics 
Statistical difference between ERQs and other QSos:
```erqStat.ipynb```
### Kernel Density Estimation
 Applying KDE on the parameter space for checking the density 
behavior of ERQs in comparison to the normal blue QSOs and obtaining the 
uncertainties based on Bootstrap re-sampling:
```KDEcontour2D.m```
### Mock LOF analysis 2D
2D Single Gaussian mock:```2D-mockLOF-1G.ipynb```
2D Double Gaussian mock:```2D-mockLOF-G12.ipynb```

### Mock LOF analysis 3D
3D Single Gaussian mock:```3D-mockLOF-1G.ipynb```
3D Double Gaussian mock:```3D-mockLOF-G12.ipynb```

### Boundary selection in 2D + Median Spectra in 2D
```2D-Analysis.ipynb```
```2d-Boundary-finder.ipynb```

### Boundary selection in 3D + Median Spectra in 3D
```3D-Analysis.ipynb```
```3d-Boundary-finder.ipynb```
```isosurf_binning_Boundary.m```
```P_Intersection_finder.m```


