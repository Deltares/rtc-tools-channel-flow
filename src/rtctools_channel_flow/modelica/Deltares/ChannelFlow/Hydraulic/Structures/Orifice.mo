within  Deltares.ChannelFlow.Hydraulic.Structures;

model Orifice "Orifice that only allows flow when HQDown.H < HQUp.H"
  extends Deltares.ChannelFlow.Hydraulic.Structures.DischargeControlledStructure(Q(min=0.0));
  //This block if from rtc-tools-hydraulic-structures and works correctly if the corresponding mixin is imported
  parameter Modelica.SIunits.Length dH_max = 10.0;
  parameter Modelica.SIunits.Area area = 1.0;
  parameter Real discharge_coefficient = 0.61;
end Orifice;
