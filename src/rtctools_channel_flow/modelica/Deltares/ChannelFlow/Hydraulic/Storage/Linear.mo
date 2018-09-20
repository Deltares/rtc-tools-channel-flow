within Deltares.ChannelFlow.Hydraulic.Storage;

model Linear "Storage with linear level-storage relation"
  extends Internal.PartialStorage(HQ.H(min = H_b), V_nominal = 1 * A, V(nominal = A));
  // Surface area
  parameter Modelica.SIunits.Area A;
  // Bed level
  parameter Modelica.SIunits.Position H_b;
equation
  V = A * (HQ.H - H_b);
end Linear;
