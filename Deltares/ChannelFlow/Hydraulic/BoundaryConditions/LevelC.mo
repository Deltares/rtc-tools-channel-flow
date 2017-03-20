within Deltares.ChannelFlow.Hydraulic.BoundaryConditions;

model LevelC
  extends Deltares.ChannelFlow.Hydraulic.BoundaryConditions.Level(redeclare connector HQPort = Deltares.ChannelFlow.Interfaces.HQZCPort);
  input Modelica.SIunits.Density C;
equation
  HQ.C = C;
end LevelC;
