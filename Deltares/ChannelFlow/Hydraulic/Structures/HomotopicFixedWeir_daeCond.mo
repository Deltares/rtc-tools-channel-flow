within Deltares.ChannelFlow.Hydraulic.Structures;

model HomotopicFixedWeir_daeCond "Homotopic Two-way Fixed Weir For Both Free and Submerged Flow Regimes"
  extends Deltares.ChannelFlow.Internal.HQTwoPort;
  import SmoothMax = Deltares.Functions.SmoothMax;
  import SmoothMin = Deltares.Functions.SmoothMin;
  import SmoothAbs = Deltares.Functions.SmoothAbs;
  // Crest level
  parameter Modelica.SIunits.Position h;
  // Crest width
  parameter Modelica.SIunits.Distance width;
  // Weir flow exponent
  parameter Real exponent(unit = "1") = 1.5;
  // Tunable parameter for the magnitude of the linear flow equation
  parameter Real linear_magnitude_factor = 0.5;
  // Paramenter describing the flow regime, used for the linear equation (1 for mostly free, 0 for mostly submerged)
  parameter Real linear_free_to_submerged_factor = 1.0;
  // parameter for tuning the non-linear flow over the weir
  parameter Real submerged_flow_factor = 1.0;
  // Homotopy Parameter
  parameter Real theta;
  Modelica.SIunits.VolumeFlowRate Q;
  // Delineates where free weir flow becomes submerged weir flow
  Real submerged_flow_ratio = 1.5;
protected
  Modelica.SIunits.VolumeFlowRate linear_Q_weir;
  // Non-linearized Submerged weir flow
  Modelica.SIunits.VolumeFlowRate _nonlinear_Q_submerged_weir;
  // Non-linearized Free weir flow
  Modelica.SIunits.VolumeFlowRate _nonlinear_Q_free_weir;
  // Combined non-linear component
  Modelica.SIunits.VolumeFlowRate nonlinear_Q_weir;
  // Helper variables
  Real _one_if_free_flow;
  Modelica.SIunits.Position _higher_level_above_crest;
  Modelica.SIunits.Position _lower_level_above_crest;
  Modelica.SIunits.Position _submerged_flow_test;
equation
  // Note: a continous aproximation of the max() and abs() functions are used frequently in these equations, in the forms
  // max(a,b) = sqrt(((a) - (b)) ^ 2 + Modelica.Constants.eps) / 2 + ((a) + (b)) / 2
  // abs(a) = sqrt(a^2 + Modelica.Constants.eps)

  HQUp.Q = Q;
  // Inflow equals outflow (no storage)
  HQUp.Q + HQDown.Q = 0.0;
  // Simple linear approximation equation
  linear_Q_weir = (1 - theta) * linear_magnitude_factor * (HQUp.H - HQDown.H * (1.0 - linear_free_to_submerged_factor) - h * linear_free_to_submerged_factor) * width;
  // Helper variables
  _higher_level_above_crest = SmoothMax(SmoothMax((HQDown.H - h), (HQUp.H - h)), 0);
  _lower_level_above_crest = SmoothMax(SmoothMin((HQDown.H - h), (HQUp.H - h)), 0);
  // higher_level - submerged_flow_ratio * lower_level
  _submerged_flow_test = SmoothMax((HQUp.H - h), (HQDown.H - h)) - submerged_flow_ratio * SmoothMin((HQUp.H - h), (HQDown.H - h));
  
  // Uses equation: QsubmergedWeir[HQUp.H, HQDown.H]:= width * (HQDown.H- h) * Sqrt[2.0 * g * submerged_flow_factor * (HQUp.H- HQDown.H)]
  // NOTE: modifies it to be accurate for what ever whatever side of the weir is higher
  _nonlinear_Q_submerged_weir = width * _lower_level_above_crest * sqrt(2.0 * Modelica.Constants.g_n * submerged_flow_factor * SmoothAbs(HQUp.H - HQDown.H)) * (HQUp.H - HQDown.H) / SmoothAbs(HQUp.H - HQDown.H);
  
  // Uses equation: QfreeWeir[HQUp.H]:= 2.0 / 3.0 * width * (2.0 / 3.0 * g) ^ 0.5 * (HQUp.H - h) ^ exponent
  // NOTE: modifies the equation to be accurate for what ever whatever side of the weir is higher
  _nonlinear_Q_free_weir = 2.0 / 3.0 * width * (2.0 / 3.0 * Modelica.Constants.g_n) ^ 0.5 * _higher_level_above_crest ^ exponent * (HQUp.H - HQDown.H) / SmoothAbs(HQUp.H - HQDown.H);

  // Filter determining if flow is free or submerged
  //Embeds the logic for If[HQUp - h > submerged_flow_ratio * (HQDown - h), Freeweir, submergedWeir]... into a single function, e.g:
  _one_if_free_flow = SmoothMax(_submerged_flow_test, 0.0) / SmoothAbs(_submerged_flow_test);
  // Combine both Nonlinear flows
  nonlinear_Q_weir = theta * (_nonlinear_Q_submerged_weir * (1 - _one_if_free_flow) + _nonlinear_Q_free_weir * _one_if_free_flow);
  // Combine Linear and Nonlinear Flows
  Q = linear_Q_weir + nonlinear_Q_weir;
  annotation(Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10}), graphics = {Polygon(visible = true, origin = {-0, -16.667}, fillColor = {128, 128, 128}, fillPattern = FillPattern.Solid, lineThickness = 2, points = {{0, 66.667}, {-50, -33.333}, {50, -33.333}})}), Diagram(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10})));
end HomotopicFixedWeir_daeCond;

