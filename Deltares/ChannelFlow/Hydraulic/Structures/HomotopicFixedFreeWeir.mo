within Deltares.ChannelFlow.Hydraulic.Structures;

model HomotopicFixedFreeWeir "Homotopic One-way Fixed Weir In the Free Flow Regime"
  // Warning!!! Should consider adding constriants HQCMDown.H(max=h) and possibly HQCMUp(min=h)
  extends Deltares.ChannelFlow.Internal.HQCMTwoPort;
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
  Modelica.SIunits.VolumeFlowRate Q;
protected
  Modelica.SIunits.VolumeFlowRate linear_Q_weir;
  Modelica.SIunits.VolumeFlowRate nonlinear_Q_weir;
equation
  linear_Q_weir = (1 - theta) * ((2.0 / 3.0 * width * sqrt(2.0 / 3.0 * Modelica.Constants.g_n ) * linearization_offset ^ exponent) / linearization_offset * (HQCMUp.H - h));
  nonlinear_Q_weir = theta * (2.0 / 3.0 * width * sqrt(2.0 / 3.0 * Modelica.Constants.g_n ) * (sqrt((HQCMUp.H - h) ^ 2 + Modelica.Constants.eps)/2 + (HQCMUp.H - h) / 2) ^ exponent);
  // Combine Linear and Nonlinear Equations
  HQCMUp.Q = linear_Q_weir + nonlinear_Q_weir;
  // Inflow equals outflow
  HQCMUp.Q + HQCMDown.Q = 0.0;
  HQCMUp.Q = Q;
  annotation(Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10}), graphics = {Polygon(visible = true, origin = {-0, -16.667}, fillColor = {128, 128, 128}, fillPattern = FillPattern.Solid, lineThickness = 2, points = {{0, 66.667}, {-50, -33.333}, {50, -33.333}})}), Diagram(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10})));
end HomotopicFixedFreeWeir;
