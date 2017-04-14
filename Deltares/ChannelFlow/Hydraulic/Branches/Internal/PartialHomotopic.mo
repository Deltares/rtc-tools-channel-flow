within Deltares.ChannelFlow.Hydraulic.Branches.Internal;

partial model PartialHomotopic
  import SI = Modelica.SIunits;
  extends Deltares.ChannelFlow.Internal.HQCMTwoPort;
  extends Deltares.ChannelFlow.Internal.QForcing;
  extends Deltares.ChannelFlow.Internal.QLateral;
  // Lateral inflow. A Matrix with n_QForcing, nQLateral rows and n_level_nodes columns. Each row corresponds to a QForcing, QLateral.Q and defines the distribution of that QForcing, QLateral.Q along the Branch.
  // NOTE: To preserve mass, each row should sum to 1.0
  parameter Real QForcing_map[n_QForcing, n_level_nodes] = fill(1 / n_level_nodes, n_QForcing, n_level_nodes);
  parameter Real QLateral_map[n_QLateral, n_level_nodes] = fill(1 / n_level_nodes, n_QLateral, n_level_nodes);
  // Wind stress
  input SI.Stress wind_stress_u(nominal = 1e-1) = 0.0; // Wind stress in x (u, easting) direction (= 0 radians, 0 degrees)
  input SI.Stress wind_stress_v(nominal = 1e-1) = 0.0; // Wind stress in y (v, northing) direction (= 0.5*pi radians, 90 degrees)
  // Flow
  SI.VolumeFlowRate[n_level_nodes + 1] Q;
  // Water level
  SI.Position[n_level_nodes] H;
  // Length
  parameter SI.Distance length = 1.0;
  // Rotation
  parameter Real rotation_deg = 0.0; // Rotation of branch relative to x (u, easting) in degrees
  // Nominal depth and width for linearized pressure term and wind stress term
  parameter SI.Distance nominal_depth[n_level_nodes + 1] = fill(1.0, n_level_nodes + 1);
  parameter SI.Distance nominal_width[n_level_nodes + 1] = fill(1.0, n_level_nodes + 1);
  // Water density
  parameter SI.Density density_water = 1000.0;
  // Bottom friction coefficient for linearized friction term
  parameter Internal.BottomFrictionCoefficient friction_coefficient = 0.0;
  // Discretization options
  parameter Boolean use_inertia = true;
  parameter Integer n_level_nodes = 2;
  // Homotopy parameter
  parameter Real theta;
  // Nominal flow used in linearization
  parameter SI.VolumeFlowRate Q_nominal = 1.0;
  // Minimum value of the divisor of the friction term.  This defaults to a nonzero value,
  // so that empty reaches won't immediately yield NaN errors.  
  parameter Real min_divisor = 1e-12;
  // Substance flow rates
  SI.VolumeFlowRate M[n_level_nodes + 1, HQCMUp.medium.n_substances];
  // Substance concentrations
  SI.Density C[n_level_nodes, HQCMUp.medium.n_substances](each min = 0);
  // Nominal substance concentrations used in linearization
  parameter Real C_nominal[HQCMUp.medium.n_substances] = fill(1, HQCMUp.medium.n_substances);
protected
  SI.Stress _wind_stress;
  parameter SI.Angle rotation_rad = Modelica.Constants.D2R * rotation_deg; // Conversion to rotation in radians
  parameter SI.Distance dx = length / (n_level_nodes - 1);
  SI.Area[n_level_nodes] _cross_section;
  SI.Distance[n_level_nodes] _dxq;
  SI.VolumeFlowRate[n_QLateral] _lat = QLateral.Q;
  SI.VolumeFlowRate[n_level_nodes] _QPerpendicular_distribution = transpose(QForcing_map) * QForcing .+ transpose(QLateral_map) * _lat;
equation
  // Store boundary values into array for convenience
  Q[1] = HQCMUp.Q;
  Q[n_level_nodes + 1] = -HQCMDown.Q;
  H[1] = HQCMUp.H;
  H[n_level_nodes] = HQCMDown.H;
  M[1, :] = HQCMUp.M;
  M[n_level_nodes + 1, :] = -HQCMDown.M;
  C[1, :] = HQCMUp.C;
  C[n_level_nodes, :] = HQCMDown.C;
  // Compute q-segment lengths
  _dxq[1] = dx / 2;
  _dxq[2:n_level_nodes - 1] = fill(dx, n_level_nodes - 2);
  _dxq[n_level_nodes] = dx / 2;
  // calculate wind_stress
  _wind_stress = wind_stress_u * cos(rotation_rad) + wind_stress_v * sin(rotation_rad);
  // Momentum equation
  // Note that the equation is formulated without any divisions, to make collocation more robust.
  for section in 2:n_level_nodes loop
    // Water momentum equation
    (if use_inertia then 1 else 0) * der(Q[section]) + theta * Modelica.Constants.g_n * 0.5 * (_cross_section[section] + _cross_section[section - 1]) * (H[section] - H[section - 1]) / dx + (1 - theta) * Modelica.Constants.g_n * (nominal_width[section] * nominal_depth[section]) * (H[section] - H[section - 1]) / dx - nominal_width[section] / density_water * _wind_stress + theta * (Modelica.Constants.g_n * Q[section] * abs(Q[section]))/ (min_divisor + friction_coefficient^2 * (0.5 * (_cross_section[section] + _cross_section[section - 1]))^2 / (nominal_width[section] + (H[section] + H[section-1]))) + (1 - theta) * (abs(Q_nominal) * Modelica.Constants.g_n) / (friction_coefficient^2 * (nominal_width[section] * nominal_depth[section])^2 / (nominal_depth[section] * 2 + nominal_width[section])) * Q[section] = 0;
    // Substance transport
    if Q[section] > 0 then
      M[section, :] = theta * Q[section] .* C[section - 1, :] + (1 - theta) * (Q_nominal * C_nominal + C_nominal * (Q[section] - Q_nominal) + Q_nominal * ((C[section - 1, :] + C[section, :]) / 2 - C_nominal));
    else
      M[section, :] = theta * Q[section] .* C[section, :] + (1 - theta) * (Q_nominal * C_nominal + C_nominal * (Q[section] - Q_nominal) + Q_nominal * ((C[section - 1, :] + C[section, :]) / 2 - C_nominal));
    end if;
  end for;
  // Mass balance equations
  // Mass balance equations for same height nodes result in relation between flows on connectors. We can therefore chain branch elements.
  // Note that every mass balance is over half of the element, the cross section of which varies linearly between the cross section at the boundary and the cross section in the middle.
  for node in 1:n_level_nodes loop
    // Water mass balance
    der(_cross_section[node]) = (Q[node] - Q[node + 1] + _QPerpendicular_distribution[node]) / _dxq[node];
    // Substance mass balance
    theta * der(_cross_section[node] * C[node, :]) + (1 - theta) * (nominal_width[node] * nominal_depth[node] * der(C[node, :]) + C_nominal * der(_cross_section[node])) = (M[node, :] - M[node + 1, :]) / _dxq[node];  
  end for;
  annotation(Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10}), graphics = {Rectangle(visible = true, fillColor = {0, 255, 255}, fillPattern = FillPattern.Solid, extent = {{-60, -20}, {60, 20}})}));
end PartialHomotopic;
