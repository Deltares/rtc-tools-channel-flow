within Deltares.ChannelFlow.Hydraulic.Structures;
model TurbineGenerator
  extends DischargeControlledStructure;
  Real PowerGeneration;
  Real efficiency;
  Real dH;
equation
  dH = HQUp.H - HQDown.H;
  PowerGeneration = Q * efficiency * dH * 1000;
end TurbineGenerator;
