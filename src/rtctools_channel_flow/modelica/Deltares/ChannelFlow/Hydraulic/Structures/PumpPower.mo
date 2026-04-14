within Deltares.ChannelFlow.Hydraulic.Structures;
model PumpPower
  extends DischargeControlledStructure;
  Real PowerDemand;
  Real efficiency;
  Real dH;
equation
  dH = HQUp.H - HQDown.H;
  PowerDemand = Q * 1.0/efficiency * dH * 1000;
end PumpPower;
