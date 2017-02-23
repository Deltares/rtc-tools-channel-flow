within Deltares.ChannelFlow.Hydraulic.Storage;

model LinearZ
  extends Deltares.ChannelFlow.Hydraulic.Storage.Linear(redeclare connector HQPort = Deltares.ChannelFlow.Interfaces.HQZCPort);
  parameter Real theta;
  parameter  Modelica.SIunits.Volume  V_nominal= A * (-1.5 - H_b);
  parameter  Modelica.SIunits.Density  C_nominal = 1;  
equation
//  V*der(HQ.C)+HQ.C*der(V) = HQ.Z;
  theta*(der(V*HQ.C)) + (1-theta)  * (C_nominal*der(V) + V_nominal * der(HQ.C) )
  = HQ.Z;
end LinearZ;
