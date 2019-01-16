within Deltares.ChannelFlow.Hydraulic.Branches.Internal;

partial model PartialHomotopic
  // Note:  Homotopic path stability can only be guaranteed if semi_implicit_step_size
  // is set to the optimization step size, and if use_convective_acceleration is set
  // to false.  Consult Baayen and Piovesan, A continuation approach to nonlinear
  // model predictive control of open channel systems, 2018 for details:
  // https://arxiv.org/abs/1801.06507
  import SI = Modelica.SIunits;
  extends Deltares.ChannelFlow.Internal.HQTwoPort;
  extends Deltares.ChannelFlow.Internal.QForcing;
  extends Deltares.ChannelFlow.Internal.QLateral;
  function smooth_switch = Deltares.ChannelFlow.Internal.Functions.SmoothSwitch;
  // Lateral inflow. A Matrix with n_QForcing, nQLateral rows and n_level_nodes columns. Each row corresponds to a QForcing, QLateral.Q and defines the distribution of that QForcing, QLateral.Q along the Branch.
  // NOTE: To preserve mass, each row should sum to 1.0
  parameter Real QForcing_map[n_QForcing, n_level_nodes] = fill(1.0 / n_level_nodes, n_QForcing, n_level_nodes);
  parameter Real QLateral_map[n_QLateral, n_level_nodes] = fill(1.0 / n_level_nodes, n_QLateral, n_level_nodes);
  // Wind stress
  input SI.Stress wind_stress_u(nominal = 1e-1) = 0.0; // Wind stress in x (u, easting) direction (= 0 radians, 0 degrees)
  input SI.Stress wind_stress_v(nominal = 1e-1) = 0.0; // Wind stress in y (v, northing) direction (= 0.5*pi radians, 90 degrees)
  // Flow
  SI.VolumeFlowRate[n_level_nodes + 1] Q(each nominal = abs(Q_nominal));
  // Water level
  SI.Position[n_level_nodes] H(min = cat(1, max(H_b[1], H_b[2]), max(H_b[1:n_level_nodes - 2], max(H_b[2:n_level_nodes - 1], H_b[3:n_level_nodes])), max(H_b[n_level_nodes - 1], H_b[n_level_nodes])));
  // Array of Bottom Levels
  parameter SI.Position[n_level_nodes] H_b;
  // Length
  parameter SI.Distance length = 1.0;
  // Rotation
  parameter Real rotation_deg = 0.0; // Rotation of branch relative to x (u, easting) in degrees
  // Nominal depth and width for linearized pressure term and wind stress term
  parameter SI.Distance nominal_depth[n_level_nodes + 1] = fill(1.0, n_level_nodes + 1);
  parameter SI.Distance nominal_width[n_level_nodes + 1] = fill(1.0, n_level_nodes + 1);
  // Water density
  parameter SI.Density density_water = 1000.0;
  // Chézy or Manning bottom friction coefficient
  parameter Real friction_coefficient = 0.0;
  // Interpret 'friction_coefficient' as Manning instead of the default Chézy
  parameter Boolean use_manning = false;
  // Discretization options
  parameter Boolean use_inertia = true;
  parameter Boolean use_convective_acceleration = false;
  parameter Boolean use_upwind = true;
  parameter Integer n_level_nodes = 2;
  // Homotopy parameter
  parameter Real theta;
  // Nominal flow used in linearization
  parameter SI.VolumeFlowRate Q_nominal = 1.0;
  // Minimum value of the divisor of the friction term.  This defaults to a nonzero value,
  // so that empty reaches won't immediately yield NaN errors.
  parameter Real min_divisor = Deltares.Constants.eps;
  // Minimum value of the sabs(Q) part of the friction term.  This defaults to a nonzero value,
  // so that sabs(Q) = sqrt(Q^2 + min_abs_Q^2) is continuously differentiable for all Q.
  parameter SI.VolumeFlowRate min_abs_Q = Deltares.Constants.eps;
  // Time step size used to create a semi-implicit discretization of the friction term.
  // Zero by default, which means that a fully implicit discretization is used.
  input SI.Duration semi_implicit_step_size = 0.0;
  // Substance flow rates
  SI.VolumeFlowRate M[n_level_nodes + 1, HQUp.medium.n_substances](each nominal = 10);
  // Substance concentrations
  SI.Density C[n_level_nodes, HQUp.medium.n_substances](each min = 0, each nominal = 1);
  // Nominal substance concentrations used in linearization
  parameter Real C_nominal[HQUp.medium.n_substances] = fill(1e-3, HQUp.medium.n_substances);
protected
  SI.Stress _wind_stress;
  Real[n_level_nodes] _dQ_sq_div_Adx(each unit = "m^3/s^2");
  parameter SI.Angle rotation_rad = Deltares.Constants.D2R * rotation_deg; // Conversion to rotation in radians
  parameter SI.Distance dx = length / (n_level_nodes - 1);
  SI.Area[n_level_nodes] _friction;
  SI.Area[n_level_nodes] _cross_section;
  SI.Area[n_level_nodes + 1] _cross_sectionq;
  SI.Distance[n_level_nodes] _wetted_perimeter;
  parameter SI.Distance[n_level_nodes] _dxq = cat(1, {dx / 2}, fill(dx, n_level_nodes - 2), {dx / 2});
  SI.VolumeFlowRate[n_QLateral] _lat = QLateral.Q;
  SI.VolumeFlowRate[n_level_nodes] _QPerpendicular_distribution = transpose(QForcing_map) * QForcing .+ transpose(QLateral_map) * _lat;
equation
  // Store boundary values into array for convenience
  Q[1] = HQUp.Q;
  Q[n_level_nodes + 1] = -HQDown.Q;
  H[1] = HQUp.H;
  H[n_level_nodes] = HQDown.H;
  M[1, :] = HQUp.M;
  M[n_level_nodes + 1, :] = -HQDown.M;
  C[1, :] = HQUp.C;
  C[n_level_nodes, :] = HQDown.C;
  // Calculate wind stress
  _wind_stress = wind_stress_u * cos(rotation_rad) + wind_stress_v * sin(rotation_rad);
  // Compute cross sections for q nodes
  _cross_sectionq[1] = _cross_section[1];
  _cross_sectionq[2:n_level_nodes] = 0.5 * (_cross_section[1:n_level_nodes - 1] .+ _cross_section[2:n_level_nodes]);
  _cross_sectionq[n_level_nodes + 1] = _cross_section[n_level_nodes];
  // Compute dQ/dx for advection term
  _dQ_sq_div_Adx[1] = 0;
  for section in 2:n_level_nodes loop
    if (use_inertia and use_convective_acceleration) then
      if use_upwind then
        _dQ_sq_div_Adx[section] = theta * (
            smooth_switch(Q[section])
            * (
                Q[section] ^ 2 / (min_divisor + _cross_sectionq[section])
                - Q[section - 1] ^ 2 / (min_divisor + _cross_sectionq[section - 1])
            )
            / _dxq[section - 1]
            + (1 - smooth_switch(Q[section]))
            * (
                Q[section + 1] ^ 2 / (min_divisor + _cross_sectionq[section + 1])
                - Q[section] ^ 2 / (min_divisor + _cross_sectionq[section])
            )
            / _dxq[section]
        ) + (1 - theta) * (
            (
                Q[section + 1]
                * Q_nominal
                / (nominal_width[section + 1] * nominal_depth[section + 1])
                - Q[section - 1]
                * Q_nominal
                / (nominal_width[section - 1] * nominal_depth[section - 1])
            )
            / (_dxq[section - 1] + _dxq[section])
        );
      else
        _dQ_sq_div_Adx[section] = theta * (
            (
                Q[section + 1] ^ 2 / (min_divisor + _cross_sectionq[section + 1])
                - Q[section - 1] ^ 2 / (min_divisor + _cross_sectionq[section - 1])
            )
            / (_dxq[section - 1] + _dxq[section])
        ) + (1 - theta) * (
            (
                Q[section + 1]
                * Q_nominal
                / (nominal_width[section + 1] * nominal_depth[section + 1])
                - Q[section - 1]
                * Q_nominal
                / (nominal_width[section - 1] * nominal_depth[section - 1])
            )
            / (_dxq[section - 1] + _dxq[section])
        );
      end if;
    else
      _dQ_sq_div_Adx[section] = 0.0;
    end if;
  end for;
  // Momentum equation
  // Note that the equation is formulated without any divisions, to make collocation more robust.
  _friction[1] = 0.0;
  for section in 2:n_level_nodes loop
    if use_manning then
      _friction[section] = friction_coefficient^2 * (
        theta * (
          Q[section] * sqrt(Q[section]^2 + min_abs_Q^2) * (min_divisor + 0.5 * (_wetted_perimeter[section] - semi_implicit_step_size * der(_wetted_perimeter[section]) + _wetted_perimeter[section - 1] - semi_implicit_step_size * der(_wetted_perimeter[section - 1])))^(8.0 / 6.0) / (min_divisor + 0.5 * (_cross_section[section] - semi_implicit_step_size * der(_cross_section[section]) + _cross_section[section - 1] - semi_implicit_step_size * der(_cross_section[section - 1])))^(14.0 / 6.0)
        ) + (1 - theta) * (
          Q[section] * sqrt(Q_nominal^2 + min_abs_Q^2) * (nominal_depth[section] * 2 + nominal_width[section])^(8.0 / 6.0) / (nominal_width[section] * nominal_depth[section])^(14.0 / 6.0)
      ));
    else
      _friction[section] =
        theta * (
          Q[section] * sqrt(Q[section]^2 + min_abs_Q^2) * (0.5 * (_wetted_perimeter[section] - semi_implicit_step_size * der(_wetted_perimeter[section]) + _wetted_perimeter[section - 1] - semi_implicit_step_size * der(_wetted_perimeter[section - 1]))) / (min_divisor + friction_coefficient^2 * (0.5 * (_cross_section[section] - semi_implicit_step_size * der(_cross_section[section]) + _cross_section[section - 1] - semi_implicit_step_size * der(_cross_section[section - 1])))^2)
        ) + (1 - theta) * (
          Q[section] * sqrt(Q_nominal^2 + min_abs_Q^2) * (nominal_depth[section] * 2 + nominal_width[section]) / (friction_coefficient^2 * (nominal_width[section] * nominal_depth[section])^2)
        );
    end if;
    // Water momentum equation
    (if use_inertia then 1 else 0) * (der(Q[section]) + _dQ_sq_div_Adx[section]) + theta * Deltares.Constants.g_n * 0.5 * (_cross_section[section] + _cross_section[section - 1]) * (H[section] - H[section - 1]) / dx + (1 - theta) * Deltares.Constants.g_n * (nominal_width[section] * nominal_depth[section]) * (H[section] - H[section - 1]) / dx - nominal_width[section] / density_water * _wind_stress + Deltares.Constants.g_n * _friction[section] = 0;
    // Substance transport
    M[section, :] = theta * (smooth_switch(Q[section]) * (Q[section] .* C[section - 1, :]) + (1 - smooth_switch(Q[section])) * (Q[section] .* C[section, :])) + (1 - theta) * (Q_nominal * C_nominal + C_nominal * (Q[section] - Q_nominal) + Q_nominal * ((if Q_nominal > 0 then C[section - 1, :] else C[section, :]) - C_nominal));
  end for;
  // Mass balance equations
  // Mass balance equations for same height nodes result in relation between flows on connectors. We can therefore chain branch elements.
  // Note that every mass balance is over half of the element, the cross section of which varies linearly between the cross section at the boundary and the cross section in the middle.
  for node in 1:n_level_nodes loop
    // Water mass balance
    theta * der(_cross_section[node]) + (1 - theta) * 0.5 * (nominal_width[node + 1] + nominal_width[node]) * der(H[node]) = (Q[node] - Q[node + 1] + _QPerpendicular_distribution[node]) / _dxq[node];
    // Substance mass balance
    theta * der(_cross_section[node] * C[node, :]) + (1 - theta) * 0.5 * (nominal_width[node + 1] * nominal_depth[node + 1] + nominal_width[node] * nominal_depth[node]) * der(C[node, :]) = (M[node, :] - M[node + 1, :]) / _dxq[node];
  end for;
  annotation(Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10}), graphics = {Rectangle(visible = true, fillColor = {0, 255, 255}, fillPattern = FillPattern.Solid, extent = {{-60, -20}, {60, 20}})}));
end PartialHomotopic;
