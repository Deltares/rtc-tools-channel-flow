within Deltares.ChannelFlow.Hydraulic.BoundaryConditions;

model LevelC
  extends Deltares.ChannelFlow.Hydraulic.BoundaryConditions.Level(redeclare connector HQPort1 = Deltares.ChannelFlow.Interfaces.HQZCPort);
  input   Modelica.SIunits.Concentration  C;
equation
  HQ.C = C;
end LevelC;