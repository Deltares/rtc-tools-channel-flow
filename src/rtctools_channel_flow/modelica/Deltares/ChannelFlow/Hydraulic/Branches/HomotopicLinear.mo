within Deltares.ChannelFlow.Hydraulic.Branches;

model HomotopicLinear
  import SI = Modelica.SIunits;
  extends Internal.PartialHomotopic(nominal_depth = fill(uniform_nominal_depth, n_level_nodes + 1), nominal_width = linspace(width_up, width_down, n_level_nodes + 1), H_b = linspace(H_b_up, H_b_down, n_level_nodes));
  // Nominal depth
  parameter SI.Distance uniform_nominal_depth;
  // Upstream Width (same 'Up' as HQUp)
  parameter SI.Distance width_up; 
  // Downstream Width (same 'Down' as HQDown)
  parameter SI.Distance width_down;
  // Array of Widths
  parameter SI.Distance width[n_level_nodes] = linspace(width_up, width_down, n_level_nodes);
  // Upstream Bottom Level (same 'Up' as HQUp)
  parameter SI.Position H_b_up; 
  // Downstream Bottom Level (same 'Down' as HQDown)
  parameter SI.Position H_b_down;
equation
  // Compute cross sections
  _cross_section = width .* (H .- H_b);
  // Compute Wetted Perimeter
  _wetted_perimeter = width .+ 2.0 * (H .- H_b);
end HomotopicLinear;
