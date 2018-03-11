within Deltares.ChannelFlow.SimpleRouting.Branches;

block Muskingum "Muskingum routing"
  extends Internal.PartialMuskingum;
  parameter Modelica.SIunits.Time K = 1.E4 "Storage constant";
  parameter Internal.MuskingumWeightingFactor x = 0.2 "Weighting factor";
equation
  K_internal = K;
  x_internal = x;
end Muskingum;