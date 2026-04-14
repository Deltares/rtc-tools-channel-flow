within Deltares.ChannelFlow.Hydraulic.Structures;

model PumpUserEquationdH "PumpUserEquationdH"
  extends DischargeControlledStructure;
  Modelica.Units.SI.Position HW;
  // Note: This block introduces a new state (HW).
  // This must be set via an equation. The choice of equation is up to the user.
  // Note: This block can be extended to support other head options
  parameter Integer head_option = 2;
  // The default head option of 2 aligns with the notation used in hydraulic structures
  Modelica.Units.SI.Distance dH;
 
equation
  if head_option == 2 then
    dH = HQDown.H - HW;
  else
    dH = HQDown.H - HQUp.H;
  end if;
end PumpUserEquationdH;