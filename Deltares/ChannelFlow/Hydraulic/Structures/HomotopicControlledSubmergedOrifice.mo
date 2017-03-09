within Deltares.ChannelFlow.Hydraulic.Structures;

model HomotopicControlledSubmergedOrifice "Rectangular Orifice element for one-way, optimized submerged orifice flow"
  // This element is valid for h_up > h_b + d_max and for all h_down
  // NOTE: In the python script, add a constraint so that Q = 0 where h_up < h_b + d_max (which requires a second boolean)
  // This is a usually a valid assumption when the set point of h_up is greater than h_b + d_max
  extends Pump;
  import ContMax = Deltares.Functions.ContMax;
  // Orifice Width
  parameter Modelica.SIunits.Distance width;
  // Orifice bottom elevation
  parameter Modelica.SIunits.Distance h_b;
  // Orifice max gate opening
  parameter Modelica.SIunits.Distance d_max;
  // Orifice constant
  parameter Real C = 1.0;
  // Homotopy Parameter
  parameter Real theta;
  // Optimized Flow Rate
  input Modelica.SIunits.VolumeFlowRate q_opt(min=0.0);
  // Average Flow Rate
  parameter Modelica.SIunits.VolumeFlowRate q_avg;
protected
  // Helper Variable
  Modelica.SIunits.VolumeFlowRate q_diff(max = 0.0);
equation
  // Calculate a physical upper bound for flow thropugh the orifice. Uses the equation:
  // Q(HUp, HDown, d) = width * C * d * (2 * g * (max(HUp - HDown, 0))) ^ 0.5
  q_diff = q_opt - theta * (width * C * d_max * (2 * Modelica.Constants.g_n * ContMax(HQUp.H - HQDown.H, 0)) ^ 0.5); 
  Q = q_avg * (1.0 - theta) + q_opt;
  annotation(Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10}), graphics = {Polygon(visible = true, origin = {0, -16.667}, fillColor = {0, 128, 0}, fillPattern = FillPattern.Solid, lineThickness = 2, points = {{0, 66.667}, {-50, -33.333}, {50, -33.333}})}), Diagram(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10})));
end HomotopicControlledSubmergedOrifice;