within Deltares.ChannelFlow.SimpleRouting.Branches.Internal;

partial block PartialMuskingum
  extends Deltares.ChannelFlow.Internal.QSISO;
protected
  Modelica.SIunits.Time K_internal "Storage constant";
  Internal.MuskingumWeightingFactor x_internal "Weighting factor";
  // We don't introduce a storage state, as this would require the user to specify
  // its initial value.  We prefer to let the user specify the initial values for the
  // flows
equation
  der(K_internal * (x_internal * QIn.Q + (1 - x_internal) * QOut.Q)) = QIn.Q - QOut.Q;
  annotation(Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10}), graphics = {Polygon(visible = true, origin = {-5, -37.5}, fillColor = {128, 128, 128}, fillPattern = FillPattern.Solid, points = {{-45, -12.5}, {-45, 17.5}, {45, 7.5}, {45, -12.5}}), Line(visible = true, origin = {-40, 10}, points = {{0, 30}, {0, -30}}), Line(visible = true, origin = {30, 5.791}, points = {{0, 34.209}, {0, -34.209}}), Line(visible = true, origin = {-5, 20}, points = {{-35, 10}, {35, -10}}), Line(visible = true, origin = {-5, 10}, points = {{-35, 0}, {35, 0}}, pattern = LinePattern.Dash)}));
end PartialMuskingum;