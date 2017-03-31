within Deltares.ChannelFlow.Hydraulic.BoundaryConditions;

model DischargeM
  extends Deltares.ChannelFlow.Hydraulic.BoundaryConditions.Discharge(redeclare connector HQPort = Deltares.ChannelFlow.Interfaces.HQCMPort);
  input Modelica.SIunits.MassFlowRate M;
equation
  HQ.M + M = 0;
end DischargeM;
