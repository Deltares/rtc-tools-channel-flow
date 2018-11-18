within Deltares.ChannelFlow.SimpleRouting.Structures;

block DischargeControlledStructure "DischargeControlledStructure"
  import SI = Modelica.SIunits;
  extends Deltares.ChannelFlow.Internal.QSISO;
  // Inputs
  input SI.VolumeFlowRate Q;
equation
  QIn.Q = QOut.Q;
  QIn.Q = Q;
  annotation(Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10}), graphics = {Polygon(visible = true, origin = {0, -16.667}, fillColor = {255, 128, 0}, fillPattern = FillPattern.Solid, lineThickness = 0.25, points = {{0, 66.667}, {-50, -33.333}, {50, -33.333}})}), Diagram(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10})));
end DischargeControlledStructure;
