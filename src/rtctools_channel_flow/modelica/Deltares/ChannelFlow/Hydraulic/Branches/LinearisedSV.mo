within Deltares.ChannelFlow.Hydraulic.Branches;

model LinearisedSV
  import SI = Modelica.Units.SI;
  extends Deltares.ChannelFlow.Internal.HQTwoPort;
  extends Deltares.ChannelFlow.Internal.QForcing;
  extends Deltares.ChannelFlow.Internal.QLateral;
  // Lateral inflow. A Matrix with n_QForcing, nQLateral rows and n_level_nodes columns. Each row corresponds to a QForcing, QLateral.Q and defines the distribution of that QForcing, QLateral.Q along the Branch.
  // NOTE: To preserve mass, each row should sum to 1.0
  parameter Real QForcing_map[n_QForcing, n_level_nodes] = fill(1.0 / n_level_nodes, n_QForcing, n_level_nodes);
  parameter Real QLateral_map[n_QLateral, n_level_nodes] = fill(1.0 / n_level_nodes, n_QLateral, n_level_nodes);
  // Wind stress
  input SI.Stress wind_stress_u(nominal = 1e-1) = 0.0; // Wind stress in x (u, easting) direction (= 0 radians, 0 degrees)
  input SI.Stress wind_stress_v(nominal = 1e-1) = 0.0; // Wind stress in y (v, northing) direction (= 0.5*pi radians, 90 degrees)
  // Flow
  SI.VolumeFlowRate[n_level_nodes + 1] Q;
  SI.VolumeFlowRate[n_level_nodes + 1] Q_relative;
  // Water Level
  SI.Position[n_level_nodes] H(min = H_b);
  SI.Position[n_level_nodes] Y_relative;
  parameter SI.Position H_nominal = 1.0;
  parameter SI.Position H_nominal_down = 1.0;
  parameter SI.Position[n_level_nodes] H_nominal_vec = linspace(H_nominal, H_nominal_down, n_level_nodes); 
  // Length
  parameter SI.Distance length = 1.0;
  // Rotation
  parameter Real rotation_deg = 0.0; // Rotation of branch relative to x (u, easting) in degrees
  // Water density
  parameter SI.Density density_water = 1000.0;
  // Manning
  parameter Real friction_coefficient = 0.0;
  // Discretization options
  parameter Integer n_level_nodes = 4;
  // Nominal flow used in linearization
  parameter SI.VolumeFlowRate Q_nominal = 1.0;

  // Width
  parameter SI.Distance width = 1.0;
  // Upstream Bottom Level (same 'Up' as HQUp)
  parameter SI.Position H_b_up; 
  // Downstream Bottom Level (same 'Down' as HQDown)
  parameter SI.Position H_b_down;
  // Array of Bottom Levels
  parameter SI.Position[n_level_nodes] H_b = linspace(H_b_up, H_b_down, n_level_nodes); 

  //input SI.Duration semi_implicit_step_size = 0.0;
  //parameter SI.VolumeFlowRate[n_level_nodes + 1] Q0;
  //SI.Position[n_level_nodes] Y;
  //parameter SI.Position Y0[n_level_nodes];

  parameter SI.Distance T0[n_level_nodes];
  parameter SI.Velocity[n_level_nodes * 2 - 1] V0;
  parameter Real[n_level_nodes + 1] Delta;
  parameter Real[n_level_nodes] Gamma;
  parameter Real[n_level_nodes + 1] C0;
  

  parameter SI.VolumeFlowRate[n_QLateral] Qlateral_nominal = fill(0.0, n_QLateral);
  parameter SI.VolumeFlowRate[n_QForcing] Qforcing_nominal = fill(0.0, n_QForcing);


protected
  SI.Stress _wind_stress;
  constant Real D2R = 3.141592653590 / 180.0;
  parameter Real g_n=9.80665;
  parameter SI.Angle rotation_rad = D2R * rotation_deg; // Conversion to rotation in radians
  parameter SI.Distance dx = length / (n_level_nodes - 1);
  SI.Area[n_level_nodes] _cross_section;
  SI.VolumeFlowRate[n_QLateral] _lat = QLateral.Q .- Qlateral_nominal;
  SI.VolumeFlowRate[n_level_nodes] _QPerpendicular_distribution = transpose(QForcing_map) * (QForcing .- Qforcing_nominal) .+ transpose(QLateral_map) * _lat;
  
  // output Real test_1;

equation
  for node in 1:(n_level_nodes+1) loop
    Q[node] = Q_relative[node] + Q_nominal;
  end for;

  H = Y_relative + H_nominal_vec + H_b; //H is level, not depth
  Q_relative[1] = HQUp.Q - Q_nominal;
  Q_relative[n_level_nodes + 1] = (-HQDown.Q) - Q_nominal;
  Y_relative[1] = HQUp.H - H_nominal - H_b_up;
  Y_relative[n_level_nodes] = HQDown.H - H_nominal_down - H_b_down;
  _wind_stress = wind_stress_u * cos(rotation_rad) + wind_stress_v * sin(rotation_rad);


  
  der(Q_relative[2]) + 2 * V0[2] * ((Q_relative[3] - Q_relative[1]) / (1.5 * dx)) + (C0[2] ^ 2 - V0[2] ^ 2) * ((T0[1] + T0[2]) / 2) * ((Y_relative[2] - Y_relative[1]) / dx) + Delta[2] * Q_relative[2] - (Gamma[1] * Y_relative[1] + Gamma[2] * Y_relative[2]) / 2 - width / density_water * _wind_stress = 0;
  
  for node in 3:n_level_nodes - 1 loop
    der(Q_relative[node]) + 2 * V0[(node - 1) * 2] * ((Q_relative[node + 1] - Q_relative[node - 1]) / (2 * dx)) + (C0[node] ^ 2 - V0[(node - 1) * 2] ^ 2) * ((T0[node - 1] + T0[node]) / 2) * ((Y_relative[node] - Y_relative[node - 1]) / dx) + Delta[node] * Q_relative[node] - (Gamma[node - 1] * Y_relative[node - 1] + Gamma[node] * Y_relative[node]) / 2  - width / density_water * _wind_stress = 0;
  end for;
  
  der(Q_relative[n_level_nodes]) + 2 * V0[2 * n_level_nodes - 2] * ((Q_relative[n_level_nodes + 1] - Q_relative[n_level_nodes - 1]) / (1.5 * dx)) + (C0[n_level_nodes] ^ 2 - V0[2 * n_level_nodes - 2] ^ 2) * ((T0[n_level_nodes - 1] + T0[n_level_nodes]) / 2) * ((Y_relative[n_level_nodes] - Y_relative[n_level_nodes - 1]) / dx) + Delta[n_level_nodes] * Q_relative[n_level_nodes] - (Gamma[n_level_nodes - 1] * Y_relative[n_level_nodes - 1] + Gamma[n_level_nodes] * Y_relative[n_level_nodes]) / 2 - width / density_water * _wind_stress = 0;
  
  
  der(_cross_section[1]) + 2 * (Q_relative[2] - Q_relative[1]) / dx - 2 * _QPerpendicular_distribution[1] / dx = 0;
  for node in 2:n_level_nodes-1 loop
    // Middle heights calculated by mass conservation
    der(_cross_section[node]) + (Q_relative[node + 1] - Q_relative[node]) / dx - _QPerpendicular_distribution[node] / dx= 0;
  end for;
  // Boundary mass conservation
  der(_cross_section[n_level_nodes]) + 2 * (Q_relative[n_level_nodes + 1] - Q_relative[n_level_nodes]) / dx - 2 * _QPerpendicular_distribution[n_level_nodes] / dx = 0;


  for node in 1:n_level_nodes loop
    Y_relative[node] = _cross_section[node] / T0[node];
  end for;

  /*
  T0[1]*der(Y_relative[1]) + 2 * (Q_relative[2] - Q_relative[1]) / dx = 0;
  for node in 2:n_level_nodes-1 loop
    // Middle heights calculated by mass conservation
    T0[node]*der(Y_relative[node]) + (Q_relative[node + 1] - Q_relative[node]) / dx = 0;
  end for;
  // Boundary mass conservation
  T0[n_level_nodes]*der(Y_relative[n_level_nodes]) + 2 * (Q_relative[n_level_nodes + 1] - Q_relative[n_level_nodes]) / dx = 0;
  */

  annotation(Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10}), graphics = {Rectangle(visible = true, fillColor = {0, 255, 255}, fillPattern = FillPattern.Solid, extent = {{-60, -20}, {60, 20}})}));
end LinearisedSV;