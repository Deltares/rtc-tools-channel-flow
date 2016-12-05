within Deltares.ChannelFlow.Hydraulic.Branches;

model HomotopicLinear
  import SI = Modelica.SIunits;
  extends Internal.PartialHomotopic(nominal_depth = fill(uniform_nominal_depth, n_level_nodes + 1), nominal_width = fill(width, n_level_nodes + 1), H(each min=H_b));
  // Nominal depth
  parameter SI.Distance uniform_nominal_depth;
  // Constant Width
  parameter SI.Distance const_width = 1.0;
  // Width
  parameter SI.Distance width = fill(const_width, n_level_nodes);
  // Constant Bottom Level
  parameter SI.Position const_H_b = 0.0;
  // Bottom level
  parameter SI.Position H_b[n_level_nodes] = fill(const_H_b, n_level_nodes);
equation
  // Compute cross sections
  _cross_section = width .* (H - H_b);
  annotation(Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10})));
end HomotopicLinear;
