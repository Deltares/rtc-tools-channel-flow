within Deltares.ChannelFlow.Hydraulic.BoundaryConditions;

model DischargeC
  extends Deltares.ChannelFlow.Hydraulic.BoundaryConditions.Discharge(redeclare connector HQPort = Deltares.ChannelFlow.Interfaces.HQCMPort);
  input Modelica.SIunits.Density C;
equation
  HQ.C = C;
end DischargeC;
