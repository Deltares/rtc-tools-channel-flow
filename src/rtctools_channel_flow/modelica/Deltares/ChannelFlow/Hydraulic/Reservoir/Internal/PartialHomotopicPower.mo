within Deltares.ChannelFlow.Hydraulic.Reservoir.Internal;

partial model PartialHomotopicPower
// this partial model adds the hydropower equation to the HomotopicVolume partial model. 
// the head difference is dynamically computed as difference between reservoir stage and tailwater level. 
// the turbine efficiency is constant. 
// The optimization starts with a linear approximation of the power equation: head difference is assumed to be constant. This value must be provided as parameter for the node.  
  import SI = Modelica.Units.SI;
  extends PartialHomotopicVolume(theta = theta);
  // Parameters
  parameter Real theta = theta;
  // head difference for the linear approximation for theta = 0
  parameter Real dH_0;
  // The efficiency term is the product of density of water, the gravity acceleration constant g and turbine efficiency (equation below). This term is assumed to be constant and provided as parameter within the model via Python. 
  parameter Real efficiency_term;
  // The following variables are not used, we use the efficiency term instead. We keep them in for historic reasons and potential later extension. 
  // parameter SI.Density density_water = 1000.0;
  // parameter Real turbine_efficiency;
  // The tailwater level depends on the total outflow. Such an equation must be specified in another model that inherits this model. 
  SI.Position H_tw; // there is no equation for H_tw here. This creates an unbalanced model. 
  // head difference Delta H
  SI.Position dH;
  // Hydro power generation
  Real Power(nominal = power_nominal);
  parameter Real power_nominal;

equation
  // Delta H equation
  dH = H - H_tw;
  // computation of the efficiency_term (a constant term in the power equation). Not used, see above. 
  // efficiency_term = turbine_efficiency * Deltares.Constants.g_n * density_water;
  // Power equation.
  Power / power_nominal = ((1 - theta) * (efficiency_term * Q_turbine * dH_0) + theta * (efficiency_term * Q_turbine * dH)) / power_nominal;
end PartialHomotopicPower;
