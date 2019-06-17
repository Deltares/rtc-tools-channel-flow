within Deltares.ChannelFlow.SimpleRouting.Branches;

block Muskingum "Muskingum routing"
  extends Internal.PartialMuskingum(K_internal=K, x_internal=x);
  parameter Modelica.SIunits.Time K = 1.E4 "Storage constant";
  parameter Internal.MuskingumWeightingFactor x = 0.2 "Weighting factor";
end Muskingum;