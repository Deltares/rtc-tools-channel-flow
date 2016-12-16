within Deltares.ChannelFlow.Hydraulic.Structures;

model SubmergedOrifice "SubmergedOrifice"
  extends Pump;
  // Orifice Width
  parameter Modelica.SIunits.Distance width;
  // Orifice hight
  parameter Modelica.SIunits.Distance height;
  // Orifice constant
  parameter Real C = 1.0;
protected
  Real LHS;
equation
  LHS = ((Q / (width * height * C)) ^ 2) / (Modelica.Constants.g_n * 2) + HQDown.H - HQUp.H;
  annotation(Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10}), graphics = {Polygon(visible = true, origin = {0, -16.667}, fillColor = {255, 128, 0}, fillPattern = FillPattern.Solid, lineThickness = 2, points = {{0, 66.667}, {-50, -33.333}, {50, -33.333}})}), Diagram(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10})));
end SubmergedOrifice;