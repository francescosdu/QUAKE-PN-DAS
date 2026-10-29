# DAS-tailored Workflow

Workflow for detection and analysis of DAS data using DAS-specific phase pickers (PhaseNet-DAS, Zhu et al. 2023).
The workflow integrates:  
- **PhaseNet-DAS** for performing **phase picking**  
- **GaMMa** (Zhu et al. 2022) for performing **phase assocation**  
- **NLLoc** (Lomax et al. 2000) for determining **hypocentral location**  

Phase picking can be performed selecting the original PhaseNet-DAS model (Zhu et al. 2023), a retrained model using earthquakes recorded in Irpinia region (Scotto di Uccio & Hysa, 2026) or a model obtained from trasfer learning of the weights of the original model (Scotto di Uccio & Hysa, 2026)

## Requirements

We suggest using a python version > 3.11.

Install the Python dependencies with:

```bash
pip install -r requirements.txt
```

## Notebooks

The **detection analysis** with PhaseNet-DAS is performed in:

`01_detection_analysis.ipynb`


This notebook integrates xdas (Trabattoni et al. 2025) for reading DAS data and applies different PhaseNet-DAS models for phase picking. In this application, we employed the pretrained PhaseNet-DAS model (Zhu et al. 2023) and local models retrained for the Irpinia region (cit ) on a low-magnitude earthquake recorded by DAS within the Irpina Near-Fault Observatory (INFO). Phase picks are associated using GaMMA (Zhu et al. 2022), ruling out false declarations and retaining coherent phase arrivals. 

**Location analysis** with NLLoc is performed in:

`02_location_analysis.ipynb`


This notebook uss the arrival times obtained from the integration of PhaseNet and GaMMA and determines the hypocentral coordinates with NLLoc (Lomax et al. 2000) using an interpolated 1D velocity model for the Irpinia region. The final location and the resulting PDF is displayed against the reference hypocentral location from INFO, along with location uncertainties and quality metrics.

**NOTE**: Before proceeding with earthquake location, travel-time grids have to be computed using "generate_NLLoc_grid.sh". P and S travel-time grid requires ~90GB of storage


**Magnitude estimation** is performed in:

`03_local_magnitude.ipynb`


This notebook evaluates the local magnitude of the events converting strain rate records in displacement measurements (Trabattoni et al. 2023) and using a local magnitude scale for the Irpinia region (Bobbio et al. 2009). Local magnitude estimation and its uncertainty is compared against the value arising from the use of seismic stations from the conventional INFO network.

## Funding
This work is part of the TRANSFORM² project which aims to improve physical and digital infrastructure across Near-Fault Observatories (NFOs) in Europe.

TRANSFORM² is funded by the European Union under project number 101188365 within the HORIZON-INFRA-2024-DEV-01-01 call.
                                                                                                                            
