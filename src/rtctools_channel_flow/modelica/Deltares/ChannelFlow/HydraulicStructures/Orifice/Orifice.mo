within Deltares.ChannelFlow.HydraulicStructures.Orifice;

model Orifice "Orifice that only allows flow when HQDown.H < HQUp.H"
  extends Deltares.ChannelFlow.Hydraulic.Structures.DischargeControlledStructure(Q(min=0.0));
  parameter Modelica.Units.SI.Length dH_max = 10.0;
  parameter Modelica.Units.SI.Area area = 1.0;
  parameter Real discharge_coefficient = 0.61;
end Orifice;
