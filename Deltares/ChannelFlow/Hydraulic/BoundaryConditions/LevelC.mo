within Deltares.ChannelFlow.Hydraulic.BoundaryConditions;

model LevelC
  extends Deltares.ChannelFlow.Hydraulic.BoundaryConditions.Level(redeclare connector HQPort = Deltares.ChannelFlow.Interfaces.HQCMPort);
  input Modelica.SIunits.Density[HQ.NOS] C;
equation
  HQ.C = C;
end LevelC;
