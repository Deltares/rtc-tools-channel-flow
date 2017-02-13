within Deltares.ChannelFlow.Hydraulic.BoundaryConditions;

model DischargeZ
  extends Deltares.ChannelFlow.Hydraulic.BoundaryConditions.Discharge(redeclare connector HQPort1 = Deltares.ChannelFlow.Interfaces.HQZCPort);
  input   Modelica.SIunits.Concentration  Z;
equation
  HQ.Z = Z;
end DischargeZ;