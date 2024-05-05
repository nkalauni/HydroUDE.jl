## Model Description:
1. GR4J model (parameters:4)
2. GR4Jsnow model (Parameters:8 { 4 additional parameters for snow melt and snow frac})
3. GR4JsnowNN Model (Parameters: 4 + NN params)

## GR4J
* optimize 4 params + initial states.
* optimize 4 params starting with optimized initial states.

## GR4JsnowNN
* optimize NN params starting with optimized 4 params from GR4J.
