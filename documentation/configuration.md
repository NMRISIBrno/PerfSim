Most of the settings are held in the `config` structure, set in the main script. Here is listed a more detailed description of all the parameters which can be set, relating to the:

- [Phantom](#phantom)
- [Time sampling](#time-sampling)
- [Arterial input function](#arterial-input-function)
- [Acquisition parameters](#acquisition-parameters) 
- [K-space sampling](#k-space-sampling)
- [Calibration scans](#calibration-scans)
- [Image reconstruction](#image-reconstruction)

## Phantom
#### `config.phantom.resizeFactor`
The phantom, passed as the input, can be scaled by changing the `resizeFactor` to a value in the range (0,1). We're keeping the value at 0.8 to prevent artifacts in the later reconstructed images. The subsampling uses nearest neighbour interpolation and the phantom is centered in the image. According to the phantom scaling, the coil sensitivities are also scaled (to keep the same FOV). 

#### `config.phantom.slices`
A 3D acquisition might be simulated, considering the input phantom as a single slice of the 3D matrix. This results in the change of the sampling period of the dynamic sequence. For a 2D acquisition, just set this parameter to 1. 

## Time sampling
#### `config.acquisition.timeAxis.scanTime`
Complete scan time (including the preconstrast phase) in seconds is set as the input. The sampling period and the number of samples is calculated automatically. Also, during the simulations, we are padding the signals with zeroes to speed up the calculation of Fourier transform.


## Arterial input function
#### `config.phantom.AIF.model`
The function name of the used AIF model is defined in the `model` parameter. Two functions were implemented - `AIF_triexpG` and `AIF_noCA`. The first AIF is a model developed at the ISI CAS for small animal hemodynamics (see [this paper](https://doi.org/10.1016/j.mri.2019.05.024)). The second function generates only zero signal (no constrast agent applied) for the simulation of preconstrast scans.

The `AIF_triexpG` model is defined by a set of parameters in the `config.phantom.AIF` field: `A`, `B`, `C`, `tau1`, `tau2`, `tau3`, `beta`. These parameters shouldn't be changed, their description can be found in the mentioned [paper](https://doi.org/10.1016/j.mri.2019.05.024). 

Generally, any model or measured AIF can be supplied for the simulations. You have to ensure that the AIF has the specific sampling period, delay and number of samples, and modify the code.  

<!---
model AIF -> add the function with the equations
measured AIF -> save as a mat file with correct sampling and comment the GenerateAIF function
--->

## Acquisition parameters
#### `config.acquisition.model`
Again, the function name of the acquisition model is defined in this parameter. We have implemented the spoiled gradient echo FLASH acquisition (`FLASH2D` function), and also the inversion recovery look-locker (`IRLL` function) sequence for calibration scans.

Other parameters and their definitions are:

<p align="center">

| Parameter | Definition |
|----------|-----------|
| `config.acquisition.TR` | Repetition time [s] |
| `config.acquisition.TE`| Echo time [s] |
| `config.acquisition.FA` | Flip angle [rad] |
| `config.acquisition.delay` | Bolus arrival time [s] |
| `config.acquisition.krho` | Scaling factor of the SI equation [-] |
| `config.acquisition.SD` | Standard deviation of the noise [-] |
| `config.acquisition.r1` | Relaxivity r<sub>1</sub> [l/mmol/s] |

</p>

In the case of `TE`, multi-echo acquisition can be simulated by giving a vector of values. 


## K-space sampling
#### `config.acquisition.kSampling.method`
The type of k-space trajectory can be set in the `method`. Three sampling types were implemented - `cartesian`, `radial` and `rosettes`, with different sets of parameters to be adjusted.

The parameters of these samplings and their definitions are:
### Cartesian
| Parameter | Definition |
|----------|-----------|
|`config.acquisition.kSampling.phEncSteps`| Number of phase encoding steps (in the y-axis dim)|
|`config.acquisition.kSampling.nSamples`| Number of echo samples (in the x-axis dim)|
|`config.acquisition.kSampling.repetitions`| Number of repetitions / frames|

### Radial
| Parameter | Definition |
|----------|-----------|
|`config.acquisition.kSampling.projections`| Number of projections|
|`config.acquisition.kSampling.nSamples`| Number of echo samples |
|`config.acquisition.kSampling.angleIncrement`| Angle increment of projections (keeping at the golden angle 111.25°)|



### Rosettes
| Parameter | Definition |
|----------|-----------|
|`config.acquisition.kSampling.Nres`| Size of the reconstructed (rectangular) image in pixels |
|`config.acquisition.kSampling.rosettes`| Number of rosettes / excitations |
|`config.acquisition.kSampling.angleIncrement`| Angle increment of rosettes (keeping at the golden angle 111.25°) |
|`config.acquisition.kSampling.nSamples`| Number of samples in a rosette |
|`config.acquisition.kSampling.n1`| Number of fast oscilations |
|`config.acquisition.kSampling.n2`| Number of slow oscilations |
| `config.acquisition.rBW` | Sampling frequency of acquisition [Hz] |
| `config.acquisition.tacq0` | Delay before acquisition of the 1st sample / 1st TE |

A rosette is constructed from multiple petals. The number of petals is given as double of the fast oscillations. To achieve a symmetric rosette, the sum of the fast and slow oscillations has to be integer. The ratio of the fast and slow oscillations sets the petal width. Angle increment sets the rotation between the whole rosettes. The size of the reconstructed image is given as a pixel width of one side of the rectangular image.

## Calibration scans
#### `model`
Three types of calibration sequences were implemented - variable TR (`vTR`), variable FA (`vFA`) and inversion recovery look-locker (`IRLL`). The type is set in the `model` parameter. 

All models are using the no-contrast-agent AIF for the simulations. Based on the k-space sampling type, the number of `repetitions` / `projections` / `rosettes` might be set. The `vTR` and `vFA` option is utilizing, just as the dynamic sequence, the FLASH acquision; the `IRLL` option uses the IRLL model. Several parameters listed below are set for the calibration, the remaining parameters are kept the same as in the case of the dynamic sequence.


### vTR and vFA 
The `vTR` and `vFA` methods are defined by a vector of parameters - TR [ms] or FA [deg] - set in the `pars` variable. 

### IRLL
The `IRLL` calibration is defined by parameters:

| Parameter | Definition |
|----------|-----------|
|`config.acquisition.tau`| Time between excitations [s] |
|`config.acquisition.td`| Time between inversion and first excitation [s] |
|`config.acquisition.tr`| Relaxation time at the end of excitation train [s] |
|`config.acquisition.nProjPerInv`| Number of projections per inversion |
|`FA`| Flip angle [rad] |
|`TE`| Echo time [s] |


## Image reconstruction
The reconstruction of the image sequence (switched on in the `reco` option) is capable of radial and Cartesian data processing. It utilizes the simple Fast Fourier transform in the case of Cartesian data, and the Non-uniform Fast Fourier transform in the case of radial data. The number of projections used for image reconstruction in the case of radial acquisition is currently set to 89 based on our experience and might be changed in the main script. The NUFFT reconstruction doesn't utilize any compressed sensing techniques and isn't recommended for sparse sampling. 


The output is a single MAT file, suitable for subsequent perfusion analysis in our [PerfLab](http://perflab.cerit-sc.cz/) software. It contains the variables `data1` with precontrast scans, `data2` with dynamic sequence, and `info` with acquisition parameters. 

