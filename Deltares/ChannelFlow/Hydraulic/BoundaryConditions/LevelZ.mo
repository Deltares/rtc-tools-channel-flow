within Deltares.ChannelFlow.Hydraulic.BoundaryConditions;

model LevelZ
  extends Deltares.ChannelFlow.Hydraulic.BoundaryConditions.Level(redeclare connector HQPort = Deltares.ChannelFlow.Interfaces.HQCMPort);
  input Modelica.SIunits.MassFlowRate M;
equation
  HQ.M + M = 0;
end LevelZ;
