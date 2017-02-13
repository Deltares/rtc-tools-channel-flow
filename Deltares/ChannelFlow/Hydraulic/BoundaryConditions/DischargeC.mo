within Deltares.ChannelFlow.Hydraulic.BoundaryConditions;

model DischargeC
  extends Deltares.ChannelFlow.Hydraulic.BoundaryConditions.Discharge(redeclare connector HQPort1 = Deltares.ChannelFlow.Interfaces.HQZCPort);
  input   Modelica.SIunits.Concentration  C;
equation
//  HQ.Z = theta * C * HQ.Q + (1-theta)*C*1.0;
  HQ.Z  + C*Q = 0;
end DischargeC;