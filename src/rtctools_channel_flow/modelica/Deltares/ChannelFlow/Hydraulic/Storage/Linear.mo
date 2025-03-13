within Deltares.ChannelFlow.Hydraulic.Storage;

model Linear "Storage with linear level-storage relation"
  /*
  Note: The default medium is FreshWater.
  To use a different medium, decalre the choice in your model file, for example
  replaceable package MyMedium = Deltares.ChannelFlow.Media.SalineWater;
  Pass this as an argument to the Linear block (redeclare package medium=MyMedium)
  */
  extends Internal.PartialStorage(HQ.H(min = H_b), V_nominal = 1 * A, V(nominal = A));
  // Surface area
  parameter Modelica.Units.SI.Area A;
  // Bed level
  parameter Modelica.Units.SI.Position H_b;
equation
  V / V_nominal = A * (HQ.H - H_b) / V_nominal;
end Linear;
