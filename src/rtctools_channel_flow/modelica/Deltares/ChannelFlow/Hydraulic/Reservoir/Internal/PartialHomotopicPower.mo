within Deltares.ChannelFlow.Hydraulic.Reservoir.Internal;

partial model PartialHomotopicPower
  import SI = Modelica.Units.SI;
  extends PartialHomotopicVolume(theta = theta);
  // Parameters
  parameter Real theta = theta;
  // Water density
  // parameter SI.Density density_water = 1000.0;
  // Parameters for power equation
  parameter Real dH_0;
  // parameter Real turbine_efficiency;
  parameter Real efficiency_term;
  // Tailwater level
  SI.Position H_tw; // no equation for H_tw, creates unbalanced model, this equation must be specified in another model that inherits this model
  // Delta H
  SI.Position dH;
  // Power
  Real Power(nominal = power_nominal);
  parameter Real power_nominal;

equation
  // Delta H equation
  dH = H - H_tw;
  // Power equation
  // efficiency_term = turbine_efficiency * Deltares.Constants.g_n * density_water;
  Power / power_nominal = ((1 - theta) * (efficiency_term * Q_turbine * dH_0) + theta * (efficiency_term * Q_turbine * dH)) / power_nominal;

end PartialHomotopicPower;
