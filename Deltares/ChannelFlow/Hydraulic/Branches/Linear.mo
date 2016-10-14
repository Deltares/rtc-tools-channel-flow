within Deltares.ChannelFlow.Hydraulic.Branches;

model Linear
  import SI = Modelica.SIunits;
  extends Internal.PartialBranch(nominal_depth = fill(uniform_nominal_depth, n_level_nodes + 1), nominal_width = fill(width, n_level_nodes + 1), H(each min=H_b));
  // Nominal depth
  parameter SI.Distance uniform_nominal_depth;
  // Width
  parameter SI.Distance width;
  // Bottom level
  parameter SI.Position H_b;
equation
  // Compute cross sections
  _cross_section = width * (H - fill(H_b, n_level_nodes));
  annotation(Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10})));
end Linear;
