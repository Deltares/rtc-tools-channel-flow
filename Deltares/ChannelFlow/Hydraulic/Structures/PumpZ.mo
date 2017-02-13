within Deltares.ChannelFlow.Hydraulic.Structures;

model PumpZ
  extends Deltares.ChannelFlow.Hydraulic.Structures.Pump(redeclare connector HQPort1 = Deltares.ChannelFlow.Interfaces.HQZCPort);
equation
  HQUp.Z + HQDown.Z = 0;
  HQUp.C = HQDown.C;
end PumpZ;