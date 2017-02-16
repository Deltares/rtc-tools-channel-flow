within Deltares.ChannelFlow.Hydraulic.Structures;

model HomotopicFixedFreeWeir "Homotopic One-way Fixed Weir In the Free Flow Regime"
  // Warning!!! Should consider adding constriants HQDown.H(max=h) and possibly HQUp(min=h)
  extends Deltares.ChannelFlow.Internal.HQTwoPort;
  // Crest level
  parameter Modelica.SIunits.Position h;
  // Crest width
  parameter Modelica.SIunits.Distance width;
  parameter Real exponent(unit = "1") = 1.5;
  // Set this to be approximately the daily maximum expected water level above the crest level
  parameter Modelica.SIunits.Distance linearization_offset = 1;
  // Homotopy Parameter
  parameter Real theta;
  // "Max Function" hardness parameter
  parameter Integer nonlinear_hardness_factor = 10;
  Modelica.SIunits.VolumeFlowRate Q;
protected
  Modelica.SIunits.VolumeFlowRate linear_Q_weir;
  Modelica.SIunits.VolumeFlowRate nonlinear_Q_weir;
equation
  linear_Q_weir = (1 - theta) * ((2.0 / 3.0 * width * sqrt(2.0 / 3.0 * Modelica.Constants.g_n) * linearization_offset ^ exponent) / linearization_offset * (HQUp.H - h));
  nonlinear_Q_weir = theta * (2.0 / 3.0 * width * sqrt(2.0 / 3.0 * Modelica.Constants.g_n) * (log(1 + exp((HQUp.H - h) * nonlinear_hardness_factor)) / nonlinear_hardness_factor) ^ exponent);
  // Combine Linear and Nonlinear Equations
  HQUp.Q = linear_Q_weir + nonlinear_Q_weir;
  // Inflow equals outflow
  HQUp.Q + HQDown.Q = 0.0;
  HQUp.Q = Q;
  annotation(Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10}), graphics = {Polygon(visible = true, origin = {-0, -16.667}, fillColor = {128, 128, 128}, fillPattern = FillPattern.Solid, lineThickness = 2, points = {{0, 66.667}, {-50, -33.333}, {50, -33.333}})}), Diagram(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10})));
end HomotopicFixedFreeWeir;