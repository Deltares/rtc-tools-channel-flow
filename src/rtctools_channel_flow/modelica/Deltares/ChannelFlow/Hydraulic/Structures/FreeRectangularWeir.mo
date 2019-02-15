within Deltares.ChannelFlow.Hydraulic.Structures;

model FreeRectangularWeir "FreeRectangularWeir"
  extends Deltares.ChannelFlow.Internal.HQTwoPort;
  // Crest level input
  input Modelica.SIunits.Position crest_level;
  // Head
  Modelica.SIunits.Distance head(min = 0.0);
  // Weir parameters
  parameter Real coefficient;
  parameter Modelica.SIunits.Distance crest_width;
  // Homotopy parameter
  parameter Real theta;
  // Minimum value of the sabs(head) part of the weir equation.  This defaults to a nonzero value,
  // so that sabs(head) = sqrt(head^2 + min_head_Q^2) is continuously differentiable for all head.
  parameter Modelica.SIunits.VolumeFlowRate min_abs_head = Deltares.Constants.eps;
  // Nominal values used in linearization
  parameter Modelica.SIunits.Distance head_nominal;
  parameter Modelica.SIunits.MassFlowRate Q_nominal = 1;
  parameter Modelica.SIunits.Density C_nominal[HQUp.medium.n_substances] = fill(1e-3, HQUp.medium.n_substances);
equation
  // Water
  head = HQUp.H - crest_level;
  HQUp.Q + HQDown.Q = 0;
  HQUp.Q = coefficient * (2.0 / 3.0) * sqrt(2.0 * Deltares.Constants.g_n) * crest_width * ((1 - theta) * (head_nominal^2 + min_abs_head^2)^0.25 + theta * (head^2 + min_abs_head^2)^0.25) * head;
  // Substances
  HQUp.M = -HQDown.M;
  // The flow is always from up to down
  HQUp.M = theta * HQUp.C * HQUp.Q + (1 - theta) * (Q_nominal * C_nominal + C_nominal * (HQUp.Q - Q_nominal) + Q_nominal * (HQUp.C - C_nominal));
  annotation(Icon(coordinateSystem( initialScale = 0.1, grid = {10, 10}), graphics = {Polygon(origin = {0, -16.67}, fillColor = {255, 128, 0}, fillPattern = FillPattern.Solid, points = {{0, 66.667}, {-50, -33.333}, {50, -33.333}, {0, 66.667}})}), Diagram(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10})));
end FreeRectangularWeir;