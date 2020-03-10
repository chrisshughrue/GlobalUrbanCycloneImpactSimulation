# Global Urban Cyclone Impact Simulation
This simulation is the basis of results in "Global spread of local cyclone damages through urban trade networks" (1). Please refer to this publication for complete description of model and results. 

## Model 
Code simulates losses of industrial output resulting from historical cyclone landfalls. Model represents direct impacts (loss owing to constraints of damaged infrastructure where cyclones make landfall) and secondary impacts (reduced output resulting from economic and material dynamics such as price spikes and supply shortages). Simulation plays out over 5 years following a direct impact, for each of ~1200 storm simulations. 

Model mechanics draw on previous studies including a supply chain damage propagation model (3) and an adaptive regional economic model (4).

## Data
Global urban industrial network (UNICORN) describes industrial exchanges among 1,686 urban areas in units of USD as a weighted, directed network. UNICORN dataset can be obtained from authors (2).  

## References

1. Shughrue, C., Werner, B.T., & Seto, K.C. (2020). Global spread of local cyclone damages through urban trade networks. Nature Sustainability. Forthcoming. 

2. [Shughrue, C., & Seto, K. C. (2018). Systemic vulnerabilities of the global urban-industrial network to hazards. Climatic change, 151(2), 173-187.](https://link.springer.com/article/10.1007/s10584-018-2293-0)

3. [Wenz, L., & Levermann, A. (2016). Enhanced economic connectivity to foster heat stress–related losses. Science advances, 2(6), e1501026.](https://advances.sciencemag.org/content/2/6/e1501026?intcmp=trendmd-adv)

4. [Hallegatte, S. (2008). An adaptive regional input‐output model and its application to the assessment of the economic cost of Katrina. Risk Analysis: An International Journal, 28(3), 779-799.](https://onlinelibrary.wiley.com/doi/abs/10.1111/j.1539-6924.2008.01046.x)

