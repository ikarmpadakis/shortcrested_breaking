# Crest height distribution 
_Model development by Al Khalili and Karmpadakis (2025) - Imperial College London_

Details are provided in Al Khalili and Karmpadakis. 2025. Breaking occurrence and dissipation in shortcrested waves in finite water depths. Coastal Engineering.  


## Crest height function - **fCrest_UAK.m**

Calculates the crest height distribution based on a transformation to second-order simulated crest heights


### Inputs
- etac_num : array of second-order simulated crest heights

- D : sea-state directionality

- Tp : sea-state peak period

- Hs : sea-state significant wave height

- D, Tp, Hs must be stated in field scale and must corrrespond to one of the cases considered by (Al Khalili and Karmpadakis, 2025)

### Outputs

- etac : array of modelled crest heights

- etac : array of modelled crest heights normalised by Hs

- Q : probability of exceedance

- bins : bins implemented in the model

- params : model parameters

## Parameters.mat 

MATLAB data file containing model coefficients:

- A : coefficients describing the amplification behaviour

- B : coefficients describing the breaking behaviour

- Mamp : mean value of ratio r of amplified waves

- Mbr : mean value of ratio r of breaking waves

- Mnbr : mean value of ratio r of non-breaking and non-amplified waves

- Pamp : probability of amplified waves

- Pbr : probability of breaking waves

- Pnbr : probability of non-breaking and non-amplified waves

- S : sea-state parameters

- bins : explanatory variable bins used in model


### The dimensions of the variables correspond to the following:

- Dimension 1 : sea-state directionality

- Dimension 2 : sea-state peak period 

- Dimension 3 : sea-state significant wave height

- Dimension 4 (if exists) : model bin as defined by variable 'bins'

