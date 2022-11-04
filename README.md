# interpolation

We aim to develop a code that performs the interpolation of spectra computed from Koester's models. Which take into account various parameters such as 
effetive temperature, surface gravity, atmospheric composition (H, He, Si..), etc.

This code employs the RegularGridInterpolator function from the scipy library, which makes the strong assumption that the grid is rectangular and regular.

