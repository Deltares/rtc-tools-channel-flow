within  Deltares.ChannelFlow.Hydraulic.Structures;

model Orifice "Orifice that only allows flow when HQDown.H < HQUp.H"
  extends Deltares.ChannelFlow.Hydraulic.Structures.DischargeControlledStructure(Q(min=0.0));
  // This block originates from rtc-tools-hydraulic-structures and requires the corresponding orifice_mixin.
  // The Orifice block represents a free-flow (non-submerged) orifice structure.
  // It can only be used in combination with the orifice_mixin, which provides
  // the necessary formulation and auxiliary variables. It uses a boolean, therefore 
  // mixed-integer optimization is needed. The flow equation is piecewise linear.
  //
  // The crest level of the orifice is fixed. Flow through the orifice is
  // constrained using a boolean variable to ensure that discharge only occurs
  // under valid hydraulic conditions (i.e. when flow is directed downhill).
  //
  // This implementation is valid only for free-flow conditions and does not
  // account for submerged behaviour.
  parameter Modelica.SIunits.Length dH_max = 10.0;
  parameter Modelica.SIunits.Area area = 1.0;
  parameter Real discharge_coefficient = 0.61;
end Orifice;
