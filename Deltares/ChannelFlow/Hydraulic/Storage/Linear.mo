within Deltares.ChannelFlow.Hydraulic.Storage;

model Linear "Storage with linear level-storage relation"
  extends Internal.PartialStorage(HQCM.H(min = H_b));
  // Surface area
  parameter Modelica.SIunits.Area A;
  // Bed level
  parameter Modelica.SIunits.Position H_b;
  parameter Real theta = 1;
  parameter Modelica.SIunits.Volume V_nominal = A * 2.0;
  parameter Modelica.SIunits.Density C_nominal[HQCM.medium.n_substances] = fill(1,HQCM.medium.n_substances);
equation
  V = A * (HQCM.H - H_b);
  // der(V*HQ.C) = HQ.M, this equation is linearized;
  theta * der(V * HQCM.C) + (1 - theta) * (C_nominal * der(V) + V_nominal * der(HQCM.C)) = HQCM.M;
end Linear;
