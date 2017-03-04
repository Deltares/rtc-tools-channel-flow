within Deltares.ChannelFlow.Hydraulic.Structures;

model BooleanOrifice "Rectangular Orifice element for both submerged and free orifice flow"
  // This element is valid for h_up > h_b + d_max and for all h_down
  // NOTE: In the python script, add a constraint so that Q = 0 where h_up < h_b + d_max (which requires a second boolean)
  // This is a usually a valid assumption when the set point of h_up is greater than h_b + d_max
  extends Pump;
  // Orifice Width
  parameter Modelica.SIunits.Distance width;
  // Orifice bottom elevation
  parameter Modelica.SIunits.Distance h_b;
  // Orifice max gate opening
  parameter Modelica.SIunits.Distance d_max;
  // Orifice constant
  parameter Real C = 1.0;
protected
  Real LHS;
equation
  // Calculate a physical upper bound for flow thropugh the orifice. Uses the equation:
  // Q(HUp, HDown, d) = width * C * d * (2 * g * (HUp - max(HDown, (h_b + d)))) ^ 0.5
  // LHS is the left-hand-side of this equation in standard form:
  // ((Q / (width * height * C)) ^ 2) / (g * 2) + max(HDown, (h_b + d)) - HQUp.H = 0
  LHS = ((Q / (width * d_max * C)) ^ 2) / (Modelica.Constants.g_n * 2) + (sqrt(((HQDown.H) - (h_b + d_max)) ^ 2 + Modelica.Constants.eps) / 2 + ((HQDown.H) + (h_b + d_max)) / 2) - HQUp.H;
  annotation(Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10}), graphics = {Polygon(visible = true, origin = {0, -16.667}, fillColor = {0, 128, 0}, fillPattern = FillPattern.Solid, lineThickness = 2, points = {{0, 66.667}, {-50, -33.333}, {50, -33.333}})}), Diagram(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10})));
end BooleanOrifice;