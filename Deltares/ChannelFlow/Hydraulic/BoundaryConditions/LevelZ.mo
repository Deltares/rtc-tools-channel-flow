within Deltares.ChannelFlow.Hydraulic.BoundaryConditions;

model LevelZ
  extends Deltares.ChannelFlow.Hydraulic.BoundaryConditions.Level(redeclare connector HQPort = Deltares.ChannelFlow.Interfaces.HQZCPort);
  input Modelica.SIunits.MassFlowRate Z;
equation
  HQ.Z + Z = 0;
end LevelZ;