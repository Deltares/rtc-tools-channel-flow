within Deltares.ChannelFlow.Hydraulic.Structures;

model HomotopicFixedWeir "Homotopic Two-way Fixed Weir For Both Free and Submerged Flow Regimes"
  extends Deltares.ChannelFlow.Internal.HQCMTwoPort;
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
  //Modelica.SIunits.VolumeFlowRate nonlinear_Q_submerged_weir;
  // Non-linearized Free weir flow
  //Modelica.SIunits.VolumeFlowRate nonlinear_Q_free_weir;
  // Non-linearized weir flow
  Modelica.SIunits.VolumeFlowRate nonlinear_Q_weir;
  // Helper variables
  //Real one_if_free_flow;
  //Modelica.SIunits.Position higher_level_above_crest;
  //Modelica.SIunits.Position lower_level_above_crest;
  //Modelica.SIunits.Position submerged_flow_test;
equation
  // Note: a continous aproximation of the max() and abs() functions are used frequently in these equations, in the forms
  // max(a,b) = sqrt(((a) - (b)) ^ 2 + Modelica.Constants.eps) / 2 + ((a) + (b)) / 2
  // abs(a) = sqrt(a^2 + Modelica.Constants.eps)

  HQCMUp.Q = Q;
  // Inflow equals outflow (no storage)
  HQCMUp.Q + HQCMDown.Q = 0.0;
  // Simple linear approximation equation
  linear_Q_weir = (1 - theta) * linear_magnitude_factor * (HQCMUp.H - HQCMDown.H * (1.0 - linear_free_to_submerged_factor) - h * linear_free_to_submerged_factor) * width;
  nonlinear_Q_weir = theta * ((width * (sqrt((-1 * (sqrt(((-(HQCMDown.H - h)) - (-(HQCMUp.H - h))) ^ 2 + Modelica.Constants.eps) / 2 + ((-(HQCMDown.H - h)) + (-(HQCMUp.H - h))) / 2)) ^ 2 + Modelica.Constants.eps) / 2 + (-1 * (sqrt(((-(HQCMDown.H - h)) - (-(HQCMUp.H - h))) ^ 2 + Modelica.Constants.eps) / 2 + ((-(HQCMDown.H - h)) + (-(HQCMUp.H - h))) / 2)) / 2) * sqrt(2.0 * Modelica.Constants.g_n * submerged_flow_factor *  sqrt((HQCMUp.H - HQCMDown.H) ^ 2 + Modelica.Constants.eps)) * (HQCMUp.H - HQCMDown.H) / sqrt((HQCMUp.H - HQCMDown.H) ^ 2 + Modelica.Constants.eps)) * (1 - ((sqrt((sqrt(((HQCMUp.H - h) - (HQCMDown.H - h)) ^ 2 + Modelica.Constants.eps) / 2 + ((HQCMUp.H - h) + (HQCMDown.H - h)) / 2 + submerged_flow_ratio * (sqrt(((-(HQCMUp.H - h)) - (-(HQCMDown.H - h))) ^ 2 + Modelica.Constants.eps) / 2 + ((-(HQCMUp.H - h)) + (-(HQCMDown.H - h))) / 2)) ^ 2 + Modelica.Constants.eps) / 2 + (sqrt(((HQCMUp.H - h) - (HQCMDown.H - h)) ^ 2 + Modelica.Constants.eps) / 2 + ((HQCMUp.H - h) + (HQCMDown.H - h)) / 2 + submerged_flow_ratio * (sqrt(((-(HQCMUp.H - h)) - (-(HQCMDown.H - h))) ^ 2 + Modelica.Constants.eps) / 2 + ((-(HQCMUp.H - h)) + (-(HQCMDown.H - h))) / 2)) / 2) / sqrt((sqrt(((HQCMUp.H - h) - (HQCMDown.H - h)) ^ 2 + Modelica.Constants.eps) / 2 + ((HQCMUp.H - h) + (HQCMDown.H - h)) / 2 + submerged_flow_ratio * (sqrt(((-(HQCMUp.H - h)) - (-(HQCMDown.H - h))) ^ 2 + Modelica.Constants.eps) / 2 + ((-(HQCMUp.H - h)) + (-(HQCMDown.H - h))) / 2)) ^ 2 + Modelica.Constants.eps))) + (2.0 / 3.0 * width * (2.0 / 3.0 * Modelica.Constants.g_n) ^ 0.5 * (sqrt((sqrt(((HQCMDown.H - h) - (HQCMUp.H - h)) ^ 2 + Modelica.Constants.eps) / 2 + ((HQCMDown.H - h) + (HQCMUp.H - h)) / 2) ^ 2 + Modelica.Constants.eps) / 2 + ((sqrt(((HQCMDown.H - h) - (HQCMUp.H - h)) ^ 2 + Modelica.Constants.eps) / 2 + ((HQCMDown.H - h) + (HQCMUp.H - h)) / 2)) / 2) ^ exponent * (HQCMUp.H - HQCMDown.H) / sqrt((HQCMUp.H - HQCMDown.H) ^ 2 + Modelica.Constants.eps)) * ((sqrt((sqrt(((HQCMUp.H - h) - (HQCMDown.H - h)) ^ 2 + Modelica.Constants.eps) / 2 + ((HQCMUp.H - h) + (HQCMDown.H - h)) / 2 + submerged_flow_ratio * (sqrt(((-(HQCMUp.H - h)) - (-(HQCMDown.H - h))) ^ 2 + Modelica.Constants.eps) / 2 + ((-(HQCMUp.H - h)) + (-(HQCMDown.H - h))) / 2)) ^ 2 + Modelica.Constants.eps) / 2 + (sqrt(((HQCMUp.H - h) - (HQCMDown.H - h)) ^ 2 + Modelica.Constants.eps) / 2 + ((HQCMUp.H - h) + (HQCMDown.H - h)) / 2 + submerged_flow_ratio * (sqrt(((-(HQCMUp.H - h)) - (-(HQCMDown.H - h))) ^ 2 + Modelica.Constants.eps) / 2 + ((-(HQCMUp.H - h)) + (-(HQCMDown.H - h))) / 2)) / 2) / sqrt((sqrt(((HQCMUp.H - h) - (HQCMDown.H - h)) ^ 2 + Modelica.Constants.eps) / 2 + ((HQCMUp.H - h) + (HQCMDown.H - h)) / 2 + submerged_flow_ratio * (sqrt(((-(HQCMUp.H - h)) - (-(HQCMDown.H - h))) ^ 2 + Modelica.Constants.eps) / 2 + ((-(HQCMUp.H - h)) + (-(HQCMDown.H - h))) / 2)) ^ 2 + Modelica.Constants.eps)));

  // Helper variables
  //Max[Max[(HQCMDown.H- h) ,(HQCMUp.H- h) ],0] 
  //higher_level_above_crest = (sqrt((sqrt(((HQCMDown.H - h) - (HQCMUp.H - h)) ^ 2 + Modelica.Constants.eps) / 2 + ((HQCMDown.H - h) + (HQCMUp.H - h)) / 2) ^ 2 + Modelica.Constants.eps) / 2 + ((sqrt(((HQCMDown.H - h) - (HQCMUp.H - h)) ^ 2 + Modelica.Constants.eps) / 2 + ((HQCMDown.H - h) + (HQCMUp.H - h)) / 2)) / 2);
  //Max[-1*Max[-(HQCMDown.H-h),-(HQCMUp.H-h)],0]
  //lower_level_above_crest = (sqrt((-1 * (sqrt(((-(HQCMDown.H - h)) - (-(HQCMUp.H - h))) ^ 2 + Modelica.Constants.eps) / 2 + ((-(HQCMDown.H - h)) + (-(HQCMUp.H - h))) / 2)) ^ 2 + Modelica.Constants.eps) / 2 + (-1 * (sqrt(((-(HQCMDown.H - h)) - (-(HQCMUp.H - h))) ^ 2 + Modelica.Constants.eps) / 2 + ((-(HQCMDown.H - h)) + (-(HQCMUp.H - h))) / 2)) / 2);
  //Max[(HQCMUp.H-h),(HQCMDown.H-h)]-submerged_flow_ratio*-1*Max[-(HQCMUp.H-h),-(HQCMDown.H-h)]
  //submerged_flow_test = (sqrt(((HQCMUp.H - h) - (HQCMDown.H - h)) ^ 2 + Modelica.Constants.eps) / 2 + ((HQCMUp.H - h) + (HQCMDown.H - h)) / 2 + submerged_flow_ratio * (sqrt(((-(HQCMUp.H - h)) - (-(HQCMDown.H - h))) ^ 2 + Modelica.Constants.eps) / 2 + ((-(HQCMUp.H - h)) + (-(HQCMDown.H - h))) / 2));

  // Uses equation: QsubmergedWeir[HQCMUp.H, HQCMDown.H]:= width * (HQCMDown.H- h) * Sqrt[2.0 * g * submerged_flow_factor * (HQCMUp.H- HQCMDown.H)]
  // NOTE: modifies it to be accurate for what ever whatever side of the weir is higher, e.g:
  // QsubmergedWeir = width * Max[-1 * Max[-(HQCMDown.H - h), -(HQCMUp.H - h)], 0] * Sqrt[2.0 * g * submerged_flow_factor * Sqrt[(HQCMUp.H - HQCMDown) ^ 2 + eps]] * (HQCMUp.H - HQCMDown.H) / Sqrt[(HQCMUp.H - HQCMDown.H) ^ 2 + eps]
  //nonlinear_Q_submerged_weir = (width * (sqrt((-1 * (sqrt(((-(HQCMDown.H - h)) - (-(HQCMUp.H - h))) ^ 2 + Modelica.Constants.eps) / 2 + ((-(HQCMDown.H - h)) + (-(HQCMUp.H - h))) / 2)) ^ 2 + Modelica.Constants.eps) / 2 + (-1 * (sqrt(((-(HQCMDown.H - h)) - (-(HQCMUp.H - h))) ^ 2 + Modelica.Constants.eps) / 2 + ((-(HQCMDown.H - h)) + (-(HQCMUp.H - h))) / 2)) / 2) * sqrt(2.0 * Modelica.Constants.g_n * submerged_flow_factor *  sqrt((HQCMUp.H - HQCMDown.H) ^ 2 + Modelica.Constants.eps)) * (HQCMUp.H - HQCMDown.H) / sqrt((HQCMUp.H - HQCMDown.H) ^ 2 + Modelica.Constants.eps));

  // Uses equation: QfreeWeir[HQCMUp.H]:= 2.0 / 3.0 * width * (2.0 / 3.0 * g) ^ 0.5 * (HQCMUp.H - h) ^ exponent,
  // NOTE: modifies the equation to be accurate for what ever whatever side of the weir is higher, e.g:
  // QfreeWeir[HQCMUp.H]:= 2.0 / 3.0 * width* (2.0 / 3.0 * g) ^0.5* Max[Max[(HQCMDown- h) ,(HQCMUp- h) ],0] ^ exponent * (HQCMUp.H- HQCMDown.H)/Sqrt[(HQCMUp.H- HQCMDown.H)^2 +eps]
  //nonlinear_Q_free_weir = (2.0 / 3.0 * width * (2.0 / 3.0 * Modelica.Constants.g_n) ^ 0.5 * (sqrt((sqrt(((HQCMDown.H - h) - (HQCMUp.H - h)) ^ 2 + Modelica.Constants.eps) / 2 + ((HQCMDown.H - h) + (HQCMUp.H - h)) / 2) ^ 2 + Modelica.Constants.eps) / 2 + ((sqrt(((HQCMDown.H - h) - (HQCMUp.H - h)) ^ 2 + Modelica.Constants.eps) / 2 + ((HQCMDown.H - h) + (HQCMUp.H - h)) / 2)) / 2) ^ exponent * (HQCMUp.H - HQCMDown.H) / sqrt((HQCMUp.H - HQCMDown.H) ^ 2 + Modelica.Constants.eps));

  // Filter determining if flow is free or submerged
  //Embeds the logic for If[HQCMUp - h > submerged_flow_ratio * (HQCMDown - h), Freeweir, submergedWeir]... into a single function, e.g:
  // Max[Max[(HQCMUp.H-h),(HQCMDown.H-h)]-submerged_flow_ratio*-1*Max[-(HQCMUp.H-h),-(HQCMDown.H-h)],0]/Abs[Max[(HQCMUp.H-h),(HQCMDown.H-h)]-submerged_flow_ratio*-1*Max[-(HQCMUp.H-h),-(HQCMDown.H-h)]]
  //one_if_free_flow = ((sqrt((sqrt(((HQCMUp.H - h) - (HQCMDown.H - h)) ^ 2 + Modelica.Constants.eps) / 2 + ((HQCMUp.H - h) + (HQCMDown.H - h)) / 2 + submerged_flow_ratio * (sqrt(((-(HQCMUp.H - h)) - (-(HQCMDown.H - h))) ^ 2 + Modelica.Constants.eps) / 2 + ((-(HQCMUp.H - h)) + (-(HQCMDown.H - h))) / 2)) ^ 2 + Modelica.Constants.eps) / 2 + (sqrt(((HQCMUp.H - h) - (HQCMDown.H - h)) ^ 2 + Modelica.Constants.eps) / 2 + ((HQCMUp.H - h) + (HQCMDown.H - h)) / 2 + submerged_flow_ratio * (sqrt(((-(HQCMUp.H - h)) - (-(HQCMDown.H - h))) ^ 2 + Modelica.Constants.eps) / 2 + ((-(HQCMUp.H - h)) + (-(HQCMDown.H - h))) / 2)) / 2) / sqrt((sqrt(((HQCMUp.H - h) - (HQCMDown.H - h)) ^ 2 + Modelica.Constants.eps) / 2 + ((HQCMUp.H - h) + (HQCMDown.H - h)) / 2 + submerged_flow_ratio * (sqrt(((-(HQCMUp.H - h)) - (-(HQCMDown.H - h))) ^ 2 + Modelica.Constants.eps) / 2 + ((-(HQCMUp.H - h)) + (-(HQCMDown.H - h))) / 2)) ^ 2 + Modelica.Constants.eps));

  // Combine Linear and Nonlinear Flows
  Q = linear_Q_weir + nonlinear_Q_weir;
  annotation(Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10}), graphics = {Polygon(visible = true, origin = {-0, -16.667}, fillColor = {128, 128, 128}, fillPattern = FillPattern.Solid, lineThickness = 2, points = {{0, 66.667}, {-50, -33.333}, {50, -33.333}})}), Diagram(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10})));
end HomotopicFixedWeir;


