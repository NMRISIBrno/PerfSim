You can supplement any input data, as long as they fulfill following requirements.

## Tissue phantom
The phantom has to be saved in a lossless image format (e.g. PNG, TIF; don't use JPG!), containing integer indices corresponding to the indices in the parametric table. The background should be set to zeroes.

The size of the image can be arbitrary in general, but has to be compatible with the coil sensitivities - we are using the size 1024x1024. Also, keep in mind that the spatial sampling of the phantom should be high enough.

## Parametric table
The perfusion parameters of the indexed tissues has to be defined in a table in CSV or XLSX format. First row has to contain parameter names, second row their units. Every following row represents the parameters of one ROI in the input phantom, with its index defined in the first column. 

The set of parameters to be defined and their units is: 

<p align="center">

| Perfusion parameters | Tissue parameters | Models | 
|----------|----|-------|
| Plasma flow F<sub>p</sub> [ml/min/ml] | Native longitudinal relaxation time T<sub>10</sub> [s] | Pharmacokinetic model PK
| Extraction fraction E  [-]                    | Native transversal relaxation time T<sub>20</sub>* [s] | Gradient correction model GCM
| Extravascular extracellular volume v<sub>e</sub> [ml/ml]   | Transversal relaxivity r<sub>2p</sub> [l/mmol/s] |
| Mean capillary transit time T<sub>c</sub>  [min]  | Transversal relaxivity r<sub>2e</sub> [l/mmol/s] |
|                       | Transversal relaxivity r<sub>2</sub>* [l/mmol/s] |

</p>

Pharmacokinetic (PK) model and gradient correction model (GCM) have to be valid function names. GCM must be defined only when using a PK model combining DCE and DSC methods (from the set of our models, only DCATH). 

We have implemented several PK models:

<p align="center">

| Abbreviation | Model name | Function name |
|----------|-----------|----------|
| ETK      | Extended Tofts-Kety | `ETK` |
| TH       | Tissue Homogeneity | `TH_trunc_FT` |
| ATH      | Adiabatic Approximation of the Tissue Homogeneity | `ATHv01` |
| DCATH    | Distributed Cappilary Adiabatic Tissue Homogeneity | `aprox_aaTH_DCATH_v5_DceDsc` |
| 2CU      | Two Compartment Uptake | `TwoCU_FT` |
| 2CX      | Two Compartment Exchange | `TwoCX` |

</p>

For their definition, see e.g. paper by [Sourbron]( https://doi.org/10.1002/nbm.2940) or [Koh](https://doi.org/10.1002/jmri.22795).


## Coil sensitivities
Coil sensitivities have to be saved as a complex matrix in a MAT format. Multiple channels are possible, stacked along the third dimension. The shape of an individual channel (1st and 2nd dimension) has to be the same as of the tissue phantom, and their orientation and geometry should be compatible. That means, the measured object e.g. should be placed in the area with the maximal coil sensitivity. 

