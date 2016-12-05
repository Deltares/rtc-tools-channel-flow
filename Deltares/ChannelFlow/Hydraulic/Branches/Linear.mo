within Deltares.ChannelFlow.Hydraulic.Branches;

model Linear
  import SI = Modelica.SIunits;
  extends Internal.PartialHomotopic(theta = 0.0, nominal_depth = fill(uniform_nominal_depth, n_level_nodes + 1), nominal_width = fill(width, n_level_nodes + 1), H(each min=H_b));
  // Nominal depth
  parameter SI.Distance uniform_nominal_depth;
  // Upstream Width (same 'Up' as HQUp)
  parameter SI.Distance width_up = 1.0; 
  // Downstream Width (same 'Down' as HQDown)
  parameter SI.Distance width_down = 1.0;
  // Array of Widths
   SI.Distance width_array[n_level_nodes];
  // Upstream Bottom Level (same 'Up' as HQUp)
  parameter SI.Distance H_b_up = 0.0; 
  // Downstream Bottom Level (same 'Down' as HQDown)
  parameter SI.Distance H_b_down = 0.0;
  // Array of Bottom Levels
   SI.Distance H_b_array[n_level_nodes];
equation
  // Compute cross sections
  for node in 1:n_level_nodes:
    width_array[node] = width_up + (H_b_down - width_up) * (node - 1) / (n_level_nodes - 1)
    H_b_array[node] = H_b_up + (H_b_down - H_b_up) * (node - 1) / (n_level_nodes - 1)
  end for;
  _cross_section = width_array .* (H - H_b_array);
  annotation(Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10})));
end Linear;
