within Deltares.ChannelFlow.Hydraulic.Structures;

model HomotopicFixedWeir "Homotopic Two-way Fixed Weir For Both Free and Submerged Flow Regimes"
  extends Deltares.ChannelFlow.Internal.HQTwoPort;
  // Crest level
  parameter Modelica.SIunits.Position h;
  // Crest width
  parameter Modelica.SIunits.Distance width;
  parameter Real exponent(unit = "1") = 1.5;
  // Set this to be approximately the daily maximum expected water level above the crest level
  parameter Real linearization_parameter = 1.0;
  // parameter for tuning the flow over the weir
  parameter Real submerged_flow_factor = 1.0;
  // Homotopy Parameter
  parameter Real theta;
  Modelica.SIunits.VolumeFlowRate Q;
  // Delineates where free weir flow becomes submerged weir flow
  Real submerged_flow_ratio = 1.5;
//protected
  Modelica.SIunits.VolumeFlowRate linear_Q_weir;
  // Non-linearized Submerged weir flow
  Modelica.SIunits.VolumeFlowRate nonlinear_Q_submerged_weir;
  // Non-linearized Free weir flow
  Modelica.SIunits.VolumeFlowRate nonlinear_Q_free_weir;
  // Filter determining if flow is free or submerged
  Real one_if_free_flow;
  // Helper variables
  Modelica.SIunits.Position higher_level_above_crest;
  Modelica.SIunits.Position lower_level_above_crest;
  Modelica.SIunits.Position submerged_flow_test;
equation
  // Note: a continous aproximation of the max() and abs() functions are used frequently in these equations, in the forms
  // max(a,b) = sqrt(((a) - (b)) ^ 2 + Modelica.Constants.eps) / 2 + ((a) + (b)) / 2
  // abs(a) = sqrt(a^2 + Modelica.Constants.eps)

  HQUp.Q = Q;
  // Inflow equals outflow (no storage)
  HQUp.Q + HQDown.Q = 0.0;
  // Simple linear approximation equation
  linear_Q_weir = (1 - theta) * linearization_parameter * (HQUp.H - HQDown.H);

  // Helper variables
  //Max[Max[(HQDown.H- h) ,(HQUp.H- h) ],0] 
  higher_level_above_crest = sqrt((sqrt(((HQDown.H - h) - (HQUp.H - h)) ^ 2 + Modelica.Constants.eps) / 2 + ((HQDown.H - h) + (HQUp.H - h)) / 2) ^ 2 + Modelica.Constants.eps) / 2 + ((sqrt(((HQDown.H - h) - (HQUp.H - h)) ^ 2 + Modelica.Constants.eps) / 2 + ((HQDown.H - h) + (HQUp.H - h)) / 2)) / 2;
  //Max[-1*Max[-(HQDown.H-h),-(HQUp.H-h)],0]
  lower_level_above_crest = sqrt((-1 * (sqrt(((-(HQDown.H - h)) - (-(HQUp.H - h))) ^ 2 + Modelica.Constants.eps) / 2 + ((-(HQDown.H - h)) + (-(HQUp.H - h))) / 2)) ^ 2 + Modelica.Constants.eps) / 2 + (-1 * (sqrt(((-(HQDown.H - h)) - (-(HQUp.H - h))) ^ 2 + Modelica.Constants.eps) / 2 + ((-(HQDown.H - h)) + (-(HQUp.H - h))) / 2)) / 2;
  //Max[(HQUp.H-h),(HQDown.H-h)]-submerged_flow_ratio*-1*Max[-(HQUp.H-h),-(HQDown.H-h)]
  submerged_flow_test = sqrt(((HQUp.H - h) - (HQDown.H - h)) ^ 2 + Modelica.Constants.eps) / 2 + ((HQUp.H - h) + (HQDown.H - h)) / 2 + submerged_flow_ratio * sqrt(((-(HQUp.H - h)) - (-(HQDown.H - h))) ^ 2 + Modelica.Constants.eps) / 2 + ((-(HQUp.H - h)) + (-(HQDown.H - h))) / 2;

  // Uses equation: QsubmergedWeir[HQUp.H, HQDown.H]:= width * (HQDown.H- h) * Sqrt[2.0 * g * submerged_flow_factor * (HQUp.H- HQDown.H)]
  // NOTE: modifies it to be accurate for what ever whatever side of the weir is higher, e.g:
  // QsubmergedWeir = width * Max[-1 * Max[-(HQDown.H - h), -(HQUp.H - h)], 0] * Sqrt[2.0 * g * submerged_flow_factor * Sqrt[(HQUp.H - HQDown) ^ 2 + eps]] * (HQUp.H - HQDown.H) / Sqrt[(HQUp.H - HQDown.H) ^ 2 + eps]
  nonlinear_Q_submerged_weir = width * lower_level_above_crest * sqrt(2.0 * Modelica.Constants.g_n * submerged_flow_factor *  sqrt((HQUp.H - HQDown.H) ^ 2 + Modelica.Constants.eps)) * (HQUp.H - HQDown.H) / sqrt((HQUp.H - HQDown.H) ^ 2 + Modelica.Constants.eps);

  // Uses equation: QfreeWeir[HQUp.H]:= 2.0 / 3.0 * width * (2.0 / 3.0 * g) ^ 0.5 * (HQUp.H - h) ^ exponent,
  // NOTE: modifies the equation to be accurate for what ever whatever side of the weir is higher, e.g:
  // QfreeWeir[HQUp.H]:= 2.0 / 3.0 * width* (2.0 / 3.0 * g) ^0.5* Max[Max[(HQDown- h) ,(HQUp- h) ],0] ^ exponent * (HQUp.H- HQDown.H)/Sqrt[(HQUp.H- HQDown.H)^2 +eps]
  nonlinear_Q_free_weir = 2.0 / 3.0 * width * (2.0 / 3.0 * Modelica.Constants.g_n) ^ 0.5 * higher_level_above_crest ^ exponent * (HQUp.H - HQDown.H) / sqrt((HQUp.H - HQDown.H) ^ 2 + Modelica.Constants.eps);

  //Embeds the logic for If[HQUp - h > submerged_flow_ratio * (HQDown - h), Freeweir, submergedWeir]... into a single function, e.g:
  // Max[Max[(HQUp.H-h),(HQDown.H-h)]-submerged_flow_ratio*-1*Max[-(HQUp.H-h),-(HQDown.H-h)],0]/Abs[Max[(HQUp.H-h),(HQDown.H-h)]-submerged_flow_ratio*-1*Max[-(HQUp.H-h),-(HQDown.H-h)]]
  one_if_free_flow = (sqrt(submerged_flow_test ^ 2 + Modelica.Constants.eps) / 2 + submerged_flow_test / 2) / sqrt(submerged_flow_test ^ 2 + Modelica.Constants.eps);

  // Combine Linear and Nonlinear Flows
  Q = linear_Q_weir + (nonlinear_Q_submerged_weir * (1 - one_if_free_flow) + nonlinear_Q_free_weir * one_if_free_flow) * theta;
  annotation(Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10}), graphics = {Polygon(visible = true, origin = {-0, -16.667}, fillColor = {128, 128, 128}, fillPattern = FillPattern.Solid, lineThickness = 2, points = {{0, 66.667}, {-50, -33.333}, {50, -33.333}})}), Diagram(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10})));
end HomotopicFixedWeir;





