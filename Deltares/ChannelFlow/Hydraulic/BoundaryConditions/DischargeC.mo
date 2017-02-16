within Deltares.ChannelFlow.Hydraulic.BoundaryConditions;

model DischargeC
  extends Deltares.ChannelFlow.Hydraulic.BoundaryConditions.Discharge(redeclare connector HQPort = Deltares.ChannelFlow.Interfaces.HQZCPort);
  input   Modelica.SIunits.Concentration  C;
equation
  HQ.C = C;
end DischargeC;