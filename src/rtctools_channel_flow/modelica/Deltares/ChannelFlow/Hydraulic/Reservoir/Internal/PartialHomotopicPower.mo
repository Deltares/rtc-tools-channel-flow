within Deltares.ChannelFlow.Hydraulic.Reservoir.Internal;

partial model PartialHomotopicPower
  import SI = Modelica.Units.SI;
  extends PartialHomotopicVolume(theta = theta);
  // Parameters
  // Water density
  parameter SI.Density density_water = 1000.0;
  // Parameters for power equation
  parameter Real dH_0;
  parameter Real turbine_efficiency;
  parameter Real power_nominal;
  parameter Real efficiency_parameter = turbine_efficiency * Deltares.Constants.g_n * density_water;
  // Tailwater level
  SI.Position H_tw;
  // Delta H
  SI.Position dH;
  // Power
  Real Power (nominal = power_nominal);
equation
  // Delta H equation
  dH = H - H_tw;
  // Power equation
  Power / power_nominal = ((1 - theta) * (efficiency_parameter * Q_turbine * dH_0) + theta * (efficiency_parameter * Q_turbine * dH)) / power_nominal;
end PartialHomotopicPower;
