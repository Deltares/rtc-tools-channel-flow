within Deltares.ChannelFlow.Hydraulic.Branches.Internal;

partial model PartialBranch
  import SI = Modelica.SIunits; 
  extends Deltares.ChannelFlow.Internal.HQTwoPort;
  extends Deltares.ChannelFlow.Internal.QForcing;
  // Lateral inflow points, use -1 for diffuse inflow
  parameter Integer QForcing_chainage[n_QForcing] = fill(1, n_QForcing);
  // Wind stress
  input SI.Stress wind_stress(nominal = 1e-1) = 0.0;
  // Flow
  SI.VolumeFlowRate[n_level_nodes + 1] Q;
  // Water level
  SI.Position[n_level_nodes] H;
  // Cross section
  SI.Area[n_level_nodes] cross_section;
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
protected
  parameter SI.Distance dx = length / (n_level_nodes - 1);
  SI.Distance[n_level_nodes] dxq;
  SI.VolumeFlowRate[n_level_nodes] QForcing_distribution;
algorithm
  // Fill helper array for lateral forcing
  QForcing_distribution := fill(0.0, n_level_nodes);
  for i_Qforcing in 1:n_QForcing loop
    if QForcing_chainage[i_Qforcing] == -1 then
        QForcing_distribution := QForcing_distribution + QForcing[i_Qforcing] * dxq / length;
    else
        QForcing_distribution[QForcing_chainage[i_Qforcing]] := QForcing_distribution[QForcing_chainage[i_Qforcing]] + QForcing[i_Qforcing];
    end if;
  end for;
equation
  // Store boundary values into array for convenience
  Q[1] = HQUp.Q;
  Q[n_level_nodes + 1] = -HQDown.Q;
  H[1] = HQUp.H;
  H[n_level_nodes] = HQDown.H;
  // Compute q-segment lengths
  dxq[1] = dx / 2;
  dxq[2:n_level_nodes - 1] = fill(dx, n_level_nodes - 2);
  dxq[n_level_nodes] = dx / 2;
  // Momentum equation
  // Note that the equation is formulated without any divisions, to make collocation more robust.
  for section in 2:n_level_nodes loop
    (if use_inertia then 1 else 0) * der(Q[section]) + Modelica.Constants.g_n * nominal_depth[section] * nominal_width[section] * (H[section] - H[section - 1]) / dx - nominal_width[section] / density_water * wind_stress + friction_coefficient * Q[section] = 0;
  end for;
  // Mass balance equations
  // Mass balance equations for same height nodes result in relation between flows on connectors.  We can therefore chain branch elements.
  // Note that every mass balance is over half of the element, the cross section of which varies linearly between the cross section at the boundary and the cross section in the middle.
  for node in 1:n_level_nodes loop
    der(cross_section[node]) = (Q[node] - Q[node + 1] + QForcing_distribution[node]) / dxq[node];
  end for;
  annotation(Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10}), graphics = {Rectangle(visible = true, fillColor = {0, 255, 255}, fillPattern = FillPattern.Solid, extent = {{-60, -20}, {60, 20}})}));
end PartialBranch;

