within Deltares.ChannelFlow.Salt.Elements;

block SaltyLinear "Storage with linear level-storage relation"
  extends Deltares.ChannelFlow.Salt.Elements.SaltyPartialStorage(HQ.H(min = H_b), V_nominal = 1 * A, V(nominal = A));
  // Surface area
  parameter Modelica.SIunits.Area A;
  // Bed level
  parameter Modelica.SIunits.Position H_b;
equation  
  V / V_nominal = A * (HQ.H - H_b) / V_nominal;
end SaltyLinear;
