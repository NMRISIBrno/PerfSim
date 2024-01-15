# Generation of ground truth data

Ground truth data for perfusion analysis can be created using the `GetParameterMaps` function (located in `matlab/Utils`). 

## Input

| Parameter | Definition | Format |
|----------|-----------|----------|
| `tissue` | parametric table (from output file `tissue.mat`) | table |
| `outputVars` | ground truth maps/signals | cell / string |
| `outputSize` | (optional) the size of output maps | integer vector (row and col size) |
| `outputTs` | (optional) the time sampling of output signals [s] | float |

If `outputSize` or `outputTs` is not given, the original size or sampling period is kept for the data. Spatial resampling is done using the nearest neighbour method, temporal resampling is done the linear interpolation.

The `outputVars` can be given as a cell with parameter names, or using one of the string commands:

- `scalar` - perfusion maps 
- `dynamic` - image sequences 
- `all` - perfusion maps and image sequences. 

All the possible variables are listed here:

<p align="center">

| Perfusion maps | | Image sequences |  |
|---|---|----|----|
| F<sub>p</sub> [ml/min/ml] | plasma flow | TRF [1/min] | tissue residue function
| E [-] | extraction fraction | c<sub>t</sub> [mmol/l] | tissue contrast agent concentration
| v<sub>e</sub> [ml/ml] | extracellular extravascular volume| R<sub>1</sub> [1/s] | longitudinal relaxation rate
| T<sub>c</sub> [min] | mean capillary transit time | R<sub>2</sub> [1/s] | transversal relaxation rate
| v<sub>p</sub> [-] | plasma volume | SI<sub>TE0</sub> [-] | signal intensity with TE=0
| K<sup>trans</sup> [1/min] | transfer constant plasma-EES| SI [-] | signal intensity
| k<sub>ep</sub> [1/min] | transfer constant EES-plasma |
| PS [ml/min/ml] | permeability-surface product |

</p>


## Output
| Parameter | Definition | Format |
|----------|-----------|----------|
| `outputData` | one-row table with selected variables and their units| table |

The units of the parameters are saved in the table properties and their listing might be accessed by the command: `outputData.Properties.VariableUnits`, similarly as all the variable names: `outputData.Properties.VariableNames`.

The sampling period of the sequence is saved as a custom property and can be accessed as `outputData.Properties.CustomProperty.SamplingPeriod`, with its unit in `outputData.Properties.CustomProperty.SamplingPeriodUnit`.
