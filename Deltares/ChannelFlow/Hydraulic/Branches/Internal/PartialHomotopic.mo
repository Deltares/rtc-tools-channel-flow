within Deltares.ChannelFlow.Hydraulic.Branches.Internal;

partial model PartialHomotopic
  import SI = Modelica.SIunits;
  extends Deltares.ChannelFlow.Internal.HQTwoPort;
  extends Deltares.ChannelFlow.Internal.QForcing;
  // Lateral inflow. A Matrix with n_QForcing rows and n_level_nodes columns. Each row corresponds to a QForcing and defines the distribution of that QForcing along the Branch.
  // NOTE: To preserve mass, each row should sum to 1.0
  parameter Real QForcing_map[n_QForcing, n_level_nodes] = fill(1 / n_level_nodes, n_QForcing, n_level_nodes);
  // Wind stress
  input SI.Stress wind_stress(nominal = 1e-1) = 0.0;
  // Flow
  SI.VolumeFlowRate[n_level_nodes + 1] Q;
  // Water level
  SI.Position[n_level_nodes] H;
  // Length
  parameter SI.Distance length = 1.0;
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
  parameter Real theta;
  parameter SI.VolumeFlowRate Q_nominal = 1.0;
protected
  parameter SI.Distance dx = length / (n_level_nodes - 1);
  SI.Area[n_level_nodes] _cross_section;
  SI.Distance[n_level_nodes] _dxq;
  SI.VolumeFlowRate[n_level_nodes] _QForcing_distribution = transpose(QForcing_map) * QForcing;
equation
  // Store boundary values into array for convenience
  Q[1] = HQUp.Q;
  Q[n_level_nodes + 1] = -HQDown.Q;
  H[1] = HQUp.H;
  H[n_level_nodes] = HQDown.H;
  // Compute q-segment lengths
  _dxq[1] = dx / 2;
  _dxq[2:n_level_nodes - 1] = fill(dx, n_level_nodes - 2);
  _dxq[n_level_nodes] = dx / 2;
  // Momentum equation
  // Note that the equation is formulated without any divisions, to make collocation more robust.
  for section in 2:n_level_nodes loop
      (if use_inertia then 1 else 0) * der(Q[section]) + theta * Modelica.Constants.g_n * 0.5 * (_cross_section[section] + _cross_section[section - 1]) * (H[section] - H[section - 1]) / dx + (1 - theta) * Modelica.Constants.g_n * (nominal_width[section] * nominal_depth[section]) * (H[section] - H[section - 1]) / dx - nominal_width[section] / density_water * wind_stress + theta * (Modelica.Constants.g_n * Q[section] * abs(Q[section]))/ (friction_coefficient^2 * (0.5 * (_cross_section[section] + _cross_section[section - 1]))^2 / (nominal_width[section] + (H[section] + H[section-1]))) + (1 - theta) * (abs(Q_nominal) * Modelica.Constants.g_n) / (friction_coefficient^2 * (nominal_width[section] * nominal_depth[section])^2 / (nominal_depth[section] * 2 + nominal_width[section])) * Q[section] = 0;
  end for;
  // Mass balance equations
  // Mass balance equations for same height nodes result in relation between flows on connectors. We can therefore chain branch elements.
  // Note that every mass balance is over half of the element, the cross section of which varies linearly between the cross section at the boundary and the cross section in the middle.
  for node in 1:n_level_nodes loop
    der(_cross_section[node]) = (Q[node] - Q[node + 1] + _QForcing_distribution[node]) / _dxq[node];
  end for;
  annotation(Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10}), graphics = {Rectangle(visible = true, fillColor = {0, 255, 255}, fillPattern = FillPattern.Solid, extent = {{-60, -20}, {60, 20}})}));
end PartialHomotopic;
