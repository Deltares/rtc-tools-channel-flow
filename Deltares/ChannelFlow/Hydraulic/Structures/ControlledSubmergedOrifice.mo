within Deltares.ChannelFlow.Hydraulic.Structures;

model ControlledSubmergedOrifice "Rectangular Orifice element for one-way, optimized submerged orifice flow"
  extends Pump;//(Q(min=0.0));
  import SmoothMax = Deltares.Functions.SmoothMax;
  // Orifice Width
  parameter Modelica.SIunits.Distance width;
  // Orifice max gate opening
  parameter Modelica.SIunits.Distance d_max;
  // Orifice constant
  parameter Real C = 1.0;
protected
  // Helper Variable
  //Modelica.SIunits.VolumeFlowRate q_diff;//(max = 0.0);
equation
  // Calculate a physical upper bound for flow thropugh the orifice. Uses the equation:
  // Q(HUp, HDown, d) = width * C * d * (2 * g * (max(HUp - HDown, 0))) ^ 0.5
  //q_diff = Q - (width * C * d_max * (2 * Modelica.Constants.g_n * SmoothMax(HQUp.H - HQDown.H, 0)) ^ 0.5); 
  annotation(Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10}), graphics = {Polygon(visible = true, origin = {0, -16.667}, fillColor = {0, 128, 0}, fillPattern = FillPattern.Solid, lineThickness = 2, points = {{0, 66.667}, {-50, -33.333}, {50, -33.333}})}), Diagram(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10})));
end ControlledSubmergedOrifice;