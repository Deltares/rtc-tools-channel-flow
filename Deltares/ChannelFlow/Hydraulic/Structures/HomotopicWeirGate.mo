within Deltares.ChannelFlow.Hydraulic.Structures;

class HomotopicWeirGate "Gate in River that follows one-way submerged weir flow"
  extends Deltares.ChannelFlow.Internal.HQCMTwoPort;
  // Crest level
  parameter Modelica.SIunits.Position bottom_level;
  // Crest width
  parameter Modelica.SIunits.Distance width;
  // parameter for tuning the linear flow over the weir
  parameter Real linear_magnitude_factor = 0.5;
  // parameter for tuning the non-linear flow over the weir
  parameter Real submerged_flow_factor = 1.0;
  // Homotopy Parameter
  parameter Real theta;
  Modelica.SIunits.VolumeFlowRate Q;
//protected
  Modelica.SIunits.VolumeFlowRate linear_Q_weir;
  // Non-linearized Submerged weir flow
  Modelica.SIunits.VolumeFlowRate oneway_nonlinear_Q_submerged_weir;
  
equation
  // Note: a continous aproximation of the max() and abs() functions are used frequently in these equations, in the forms
  // max(a,b) = sqrt(((a) - (b)) ^ 2 + Modelica.Constants.eps) / 2 + ((a) + (b)) / 2
  // abs(a) = sqrt(a^2 + Modelica.Constants.eps)

  HQCMUp.Q = Q;
  // Inflow equals outflow (no storage)
  HQCMUp.Q + HQCMDown.Q = 0.0;
  // Simple linear approximation equation
  linear_Q_weir = (1 - theta) * linear_magnitude_factor * (HQCMUp.H - HQCMDown.H) * width;
  oneway_nonlinear_Q_submerged_weir = theta * width * (HQCMDown.H - bottom_level) * sqrt(2.0 * Modelica.Constants.g_n * submerged_flow_factor *  sqrt((HQCMUp.H - HQCMDown.H) ^ 2 + Modelica.Constants.eps)) * (sqrt(((HQCMUp.H - HQCMDown.H) / sqrt((HQCMUp.H - HQCMDown.H) ^ 2 + Modelica.Constants.eps)) ^ 2 + Modelica.Constants.eps) / 2 + ((HQCMUp.H - HQCMDown.H) / sqrt((HQCMUp.H - HQCMDown.H) ^ 2 + Modelica.Constants.eps)) / 2);
  HQCMUp.Q = linear_Q_weir + oneway_nonlinear_Q_submerged_weir;
  annotation(Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10}), graphics = {Polygon(visible = true, origin = {-0, -16.667}, fillColor = {128, 128, 128}, fillPattern = FillPattern.Solid, lineThickness = 2, points = {{0, 66.667}, {-50, -33.333}, {50, -33.333}})}), Diagram(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10})));
end HomotopicWeirGate;
