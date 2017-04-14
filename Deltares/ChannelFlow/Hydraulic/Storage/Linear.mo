within Deltares.ChannelFlow.Hydraulic.Storage;

model Linear "Storage with linear level-storage relation"
  extends Internal.PartialStorage(HQCM.H(min = H_b));
  // Surface area
  parameter Modelica.SIunits.Area A;
  // Bed level
  parameter Modelica.SIunits.Position H_b;
equation
  V = A * (HQCM.H - H_b);
  // der(V*HQ.C) = HQ.M, this equation is linearized;
  theta * der(V * HQCM.C) + (1 - theta) * (C_nominal * der(V) + V_nominal * der(HQCM.C)) = HQCM.M;
end Linear;
