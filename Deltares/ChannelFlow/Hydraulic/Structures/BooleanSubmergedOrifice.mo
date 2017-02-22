within Deltares.ChannelFlow.Hydraulic.Structures;

model BooleanSubmergedOrifice "BooleanSubmergedOrifice"
  extends Pump;
  // Orifice Width
  parameter Modelica.SIunits.Distance width;
  // Orifice hight
  parameter Modelica.SIunits.Distance height;
  // Orifice constant
  parameter Real C = 1.0;
  // Overestimate of Max Flow Rate TODO: Is this needed, or can we use big M?
  parameter Modelica.SIunits.VolumeFlowRate overestimate_max_Q = 1e3;
  // Handy Big M
  parameter Real M = 1e10;
  // Boolean to detect if flow is downhill
  input Boolean is_downhill;
protected
  Real nominal_Q(min=0.0, max=1.0);
  Real is_downhill_test1(max=0.0);
  Real is_downhill_test2(min=0.0);
  Real LHS(max=0.0);
equation
  //Force Downhill Only
  nominal_Q = Q / overestimate_max_Q + (if is_downhill then 0.0 else 1.0);
  // Ensure is_downhill is true only when HQDown.H is less than HQUp.H
  is_downhill_test1 = HQDown.H - HQUp.H - (if is_downhill then 0.0 else 1.0) * M;
  is_downhill_test2 = HQDown.H - HQUp.H + (if is_downhill then 1.0 else 0.0) * M;
  // Apply a physical upper bound the flow. Uses the equation:
  // Q(HUp, HDown, d) = width * C * d * (2 * g * (HUp - HDown)) ^ 0.5
  LHS = ((Q / (width * height * C)) ^ 2) / (Modelica.Constants.g_n * 2) + HQDown.H - HQUp.H - M * (if is_downhill then 0.0 else 1.0);
  annotation(Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10}), graphics = {Polygon(visible = true, origin = {0, -16.667}, fillColor = {255, 128, 0}, fillPattern = FillPattern.Solid, lineThickness = 2, points = {{0, 66.667}, {-50, -33.333}, {50, -33.333}})}), Diagram(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10})));
end BooleanSubmergedOrifice;