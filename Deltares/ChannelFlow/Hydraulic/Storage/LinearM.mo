within Deltares.ChannelFlow.Hydraulic.Storage;

model LinearM
  extends Deltares.ChannelFlow.Hydraulic.Storage.Linear(redeclare connector HQPort = Deltares.ChannelFlow.Interfaces.HQCMPort);
  parameter Real theta;
  parameter Modelica.SIunits.Volume V_nominal = A * 2.0;
  parameter Modelica.SIunits.Density C_nominal = 1;  
equation
  //der(C*HQ.C) = HQ.Z, this equation is linearized;
  theta * der(V * HQ.C) + (1 - theta) * (C_nominal * der(V) + V_nominal * der(HQ.C)) = HQ.M;
end LinearM;
