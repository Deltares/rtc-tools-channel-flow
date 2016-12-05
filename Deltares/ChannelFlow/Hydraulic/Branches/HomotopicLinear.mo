within Deltares.ChannelFlow.Hydraulic.Branches;

model HomotopicLinear
  import SI = Modelica.SIunits;
  extends Internal.PartialHomotopic(nominal_depth = fill(uniform_nominal_depth, n_level_nodes + 1), nominal_width = fill((H_b_down + width_up) / 2, n_level_nodes + 1));
  extends Internal.InterpolatedLinearCrossSections;
  // Nominal depth
  parameter SI.Distance uniform_nominal_depth;
equation
  annotation(Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10})));
end HomotopicLinear;
