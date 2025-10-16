within Deltares.ChannelFlow.Hydraulic.Structures;
 
model Pump "Pump"
  extends DischargeControlledStructure;
  // Note: dH is a helper state which is determined via the head_option parameter
  // If dH cannot be calculated directly from the upstream and downstream heads
  // use the PumpDynamicHead model instead
  parameter Integer head_option = 0;
  Modelica.Units.SI.Distance dH;
 
equation
  if head_option == -1 then
    dH = HQUp.H;
  elseif head_option == 1 then
    dH = HQDown.H;
  else
    dH = HQDown.H - HQUp.H;
  end if;
end Pump;