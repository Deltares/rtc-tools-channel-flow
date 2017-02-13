within Deltares.ChannelFlow.Hydraulic.Storage;

model LinearZ
  extends Deltares.ChannelFlow.Hydraulic.Storage.Linear(redeclare connector HQPort1 = Deltares.ChannelFlow.Interfaces.HQZCPort);
  // Modelica.SIunits.Mass m(min=0, nominal = 1e6);
  parameter Real theta;  
equation
  V*der(HQ.C)+HQ.C*der(V) = HQ.Z;
end LinearZ;