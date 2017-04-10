within Deltares.ChannelFlow.Hydraulic.BoundaryConditions;

model LevelM
  extends Deltares.ChannelFlow.Hydraulic.BoundaryConditions.Level(redeclare connector HQPort = Deltares.ChannelFlow.Interfaces.HQCMPort);
  input Modelica.SIunits.MassFlowRate[HQ.NOS] M;
equation
  HQ.M =- M;
end LevelM;
