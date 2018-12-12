within Deltares.ChannelFlow.Hydraulic.Branches;

model HomotopicTrapezoidal
  import SI = Modelica.SIunits;
  extends Internal.PartialHomotopic(nominal_depth = fill(uniform_nominal_depth, n_level_nodes + 1), nominal_width = 0.5 * nominal_depth ./ tan(Deltares.Constants.D2R * linspace(left_slope_angle_up, left_slope_angle_down, n_level_nodes + 1)) .+ linspace(bottom_width_up, bottom_width_down, n_level_nodes + 1) .+ 0.5 * nominal_depth ./ tan(Deltares.Constants.D2R * linspace(right_slope_angle_up, right_slope_angle_down, n_level_nodes + 1)), H_b = linspace(H_b_up, H_b_down, n_level_nodes));
  // Nominal depth
  parameter SI.Distance uniform_nominal_depth;
  // Upstream Bottom Width (same 'Up' as HQUp)
  parameter SI.Distance bottom_width_up;
  // Downstream Bottom Width (same 'Down' as HQDown)
  parameter SI.Distance bottom_width_down;
  // Array of Bottom Widths
  parameter SI.Distance[n_level_nodes] bottom_width = linspace(bottom_width_up, bottom_width_down, n_level_nodes);
  // Upstream Left Slope Angle (same as 'Up' in HQUp).  Left slope = slope left when facing along the positive flow direction.
  parameter Real left_slope_angle_up(unit = "deg", min = 0.0, max = 90.0) = 90.0;
  // Downstream Left Slope Angle (same as 'Down' in HQDown).  Left slope = slope left when facing along the positive flow direction.
  parameter Real left_slope_angle_down(unit = "deg", min = 0.0, max = 90.0) = 90.0;
  // Array of Left Slope Angles.  Left slope = slope left when facing along the positive flow direction.
  parameter Real[n_level_nodes] left_slope_angle(each unit = "deg") = linspace(left_slope_angle_up, left_slope_angle_down, n_level_nodes);
  // Upstream Left Slope Angle (same as 'Up' in HQUp).  Right slope = slope right when facing along the positive flow direction.
  parameter Real right_slope_angle_up(unit = "deg", min = 0.0, max = 90.0) = 90.0;
  // Downstream Left Slope Angle (same as 'Down' in HQDown).  Right slope = slope right when facing along the positive flow direction.
  parameter Real right_slope_angle_down(unit = "deg", min = 0.0, max = 90.0) = 90.0;
  // Array of Left Slope Angles.  Right slope = slope right when facing along the positive flow direction.
  parameter Real[n_level_nodes] right_slope_angle(each unit = "deg") = linspace(right_slope_angle_up, right_slope_angle_down, n_level_nodes);
  // Upstream Bottom Level (same 'Up' as HQUp)
  parameter SI.Position H_b_up;
  // Downstream Bottom Level (same 'Down' as HQDown)
  parameter SI.Position H_b_down;
equation
  // Compute nonlinear cross sections.  These are replaced into the model by RTC-Tools (thanks to the underscore prefix),
  // into the nonlinear parts of the homotopic equations.  No separate homotopy is therefore required here.
  _cross_section = (
    0.5 * (H .- H_b) ./ tan(Deltares.Constants.D2R * left_slope_angle) .+
    bottom_width .+
    0.5 * (H .- H_b) ./ tan(Deltares.Constants.D2R * right_slope_angle)
  ) .* (H .- H_b);
  // Compute Wetted Perimeter
  _wetted_perimeter = (H .- H_b) ./ sin(Deltares.Constants.D2R .* left_slope_angle) .+ bottom_width .+ (H .- H_b) ./ sin(Deltares.Constants.D2R .* right_slope_angle);
end HomotopicTrapezoidal;
