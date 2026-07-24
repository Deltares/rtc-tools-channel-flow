HomotopicPower
==============

This model extends HomotopicVolume with a hydropower generation equation.
  

  
  
  Power generation depends on turbine discharge and the available hydraulic head,
  computed as the difference between the reservoir water level and the tailwater level.
  Turbine efficiency is assumed to be constant.
  

  
  
  The optimization starts from a simplified power equation in which the hydraulic
  head is assumed to be constant. This reference head must be provided as a model parameter.
