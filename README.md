# OroGlobo

OroGlobo (OROGraphic ancillary files generator for GLOBal atmospheric mOdels) is a collection of python scripts designed to allow the creation of the ancillary files needed by atmospheric models to force orographic drag parameterizations, namely Orographic Gravity Wave Drag (OGWD), low level flow blocking (BLOCK) and Turbulent Orographic Form Drag (TOFD) schemes. 

All orographic scales are found to influence the atmospheric flow (Kanehama et al., 2019). In fact, the  parameterization of unresolved drag have been recognized as crucial to obtain a realistic mid-latitude circulation (Sandu et al., 2019) and to reduce some of the long-standing circulation biases affecting climate
models (Pithan et al., 2016). To work properly, these schemes need suitable boundary conditions. Typically, they
need information on the physical features of unresolved orography. However, the strategies for the generation of such information, stored in models' ancillary files, can vary a lot between different modeling centers, and it is an important source of uncertainty (Elvidge et al., 2019). The procedures for the generation of subgrid orography fields is often poorly documented and software are usually not publicy available.

OroGlobo is a unique and powerful tool. It gathers the main state-of-the-art procedures for generating model mean orography and boundary conditions for ESM orographic drag parameterization in a single, open-source code, allowing the user to use state-of-the-art elevation data to build suitable ancillary files for atmospheric models. The chain of computing steps is completely transparent and follows the best practices in the field, allowing the user to have the full control over the whole procedure.

OroGlobo allows to:

- generate an operational model mean orography. The user can decide the level of additional smoothing needed, to meet numerical stability or to match the effective orographic resolution of the model dynamical core (Kanehama et al., 2019)
- generate a set of unresolved orography parameters fields, consistently representing all the scales not explicitly resolved by the operational mean orography. The parameters are those typically needed by state-of-the art parameterization schemes (e.g. Lott and Miller, 1997; Beljaars et al, 2004; Van Niekerk and Vosper, 2021; Xue and Shen, 2023):
    - unresolved orography standard deviation on scales above and below 5 km, to be used by OGWD/BLOCK and TOFD parameterizations respectively
    - unresolved orography orientation on scales above and below 5 km, to be used by OGWD/BLOCK and TOFD parameterizations respectively
    - unresolved orography slope on scales above and below 5 km, to be used by OGWD/BLOCK and TOFD parameterizations respectively
    - unresolved orography anisotropy on scales above and below 5 km, to be used by OGWD/BLOCK and TOFD parameterizations respectively
    - F1, F2, F3, hamp parameters needed by the Van Niekerk and Vosper, 2021 scheme


OroGlobo makes use of Copernicus GLO-90 elevation data [https://doi.org/10.5270/ESA-c5d3d65](https://doi.org/10.5270/ESA-c5d3d65) and produces output files in netcdf format. While it is designed to produce sub-grid orography files for the atmospheric model [GloboNe](https://git.isac.cnr.it/esm/globone) developed at [ISAC-CNR](https://www.isac.cnr.it/), it can be easily adapted to work with any regular latitude-longitude grid. At present, only this kind of grid is supported.

## Installation and requirements

The user need to create a python environment with the following libraries:

xarray (tested with version 2025.6.1)
haversine (tested with version 2.8.0)
LatLon23 (tested with version 1.0.7)
rioxarray (tested with version 0.15.0)
shapely (tested with version 2.0.2)

other than the standard python libraries numpy, scipy, matplotlib. It also requires multiprocessing to exploit parallelism and speed up calculations.

## How to use OroGlobo

The original Copernicus GLO-90 DEM data were regridded to a 1 km grid and stored in a single netcdf file. This file is needed by OroGlobo and can be downloaded [here](https://zenodo.org/records/17414359).

The code also needs:
    - information on the model grid stored in text format (examples for the GlobNe model can be found in *data_in/*)
    - information on the model land sea mask and land-water fraction in netcdf format (examples for the GlobNe model can be found in *data_in/*)

The execution of the code is divided in 5 steps, correspongind to different scripts.

- *00_oroglobo_init.py* initialize the folders sctructure
- *01_filter_1km_orog.py* filters the Copernicus 1km global orography dataset to a scale corresponding to the maximum grid spacing S of the target operational model grid. The result is a 1km resolution map, but with all scales below S filtered out. Creating this field allows to filter uniformly over the globe the orography to a resolution consistent with the maximum grid spacing S of the model target grid; in this way, in the next step we will obtain a mean orography which represents the same scales over the globe even if the target grid spacing is not uniform. This step is particularly time-consuming. GloboNe users can foun pre-filtered fields [here](https://zenodo.org/records/17435467), so they can skip this step.
- *02_from_1km_smooth_to_model_grid.oy*: takes the 1km reoslution orography smoothed to target model grid scale
    and aggregate to target model grid boxes
- *03_make_operational_orog.oy*: eventually applies additional smoothing to the model gridbox-averaged orography, and applies the model land-sea mask
- *04_make_ogwd_and_tofd_parameters.py*: calculates the orographic parameters needed by parameterizations
- *05_create_model_files.pt*: performs the last operations needed to store OroGlobo output in netcdf files with a structure fully compatible with the GloboNe model.

The file *oroglobo_parameters.yaml* must be used to configure the behaviour of OroGlobo and to set the proper paths on the user system.
The path to the 1km orographic dataset must be set through the *paths_in --> copernicus_lowresnc_global* parameter

## References



## License

Copyright 2025 Guido Davoli

OroGlobo is free software made available under the MIT License. For details see the LICENSE file.


### Fincancial Support

Financial support from ICSC – Centro Nazionale di Ricerca in High Performance Computing, Big Data and Quantum Computing, funded by European Union – NextGenerationEU

#### Copernicus Data Copyright notice:
*Produced using Copernicus WorldDEM-90 © DLR e.V. 2010-2014 and © Airbus Defence and Space GmbH 2014-2018 provided under COPERNICUS by the European Union and ESA; all rights reserved*