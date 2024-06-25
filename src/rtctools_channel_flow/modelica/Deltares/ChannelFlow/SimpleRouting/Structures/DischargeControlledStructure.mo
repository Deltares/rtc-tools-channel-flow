within Deltares.ChannelFlow.SimpleRouting.Structures;

block DischargeControlledStructure "DischargeControlledStructure"
  import SI = Modelica.Units.SI;
  extends Deltares.ChannelFlow.Internal.QSISO(QIn.Q(nominal=Q_nominal), QOut.Q(nominal=Q_nominal));
  // Inputs
  input SI.VolumeFlowRate Q;
  parameter SI.VolumeFlowRate Q_nominal = 1.0;
equation
  QIn.Q / Q_nominal = QOut.Q / Q_nominal;
  QIn.Q / Q_nominal = Q / Q_nominal;
  annotation(Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10}), graphics = {Polygon(visible = true, origin = {0, -16.667}, fillColor = {255, 128, 0}, fillPattern = FillPattern.Solid, lineThickness = 0.25, points = {{0, 66.667}, {-50, -33.333}, {50, -33.333}})}), Diagram(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10})));
end DischargeControlledStructure;
