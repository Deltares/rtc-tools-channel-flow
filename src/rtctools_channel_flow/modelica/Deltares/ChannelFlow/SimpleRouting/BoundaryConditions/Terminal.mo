within Deltares.ChannelFlow.SimpleRouting.BoundaryConditions;

block Terminal
  import SI = Modelica.Units.SI;
  extends Deltares.ChannelFlow.Internal.QSI(QIn.Q(nominal=Q_nominal));
  // Outputs
  output SI.VolumeFlowRate Q(nominal=Q_nominal);

  parameter SI.VolumeFlowRate Q_nominal = 1.0;
equation
  Q / Q_nominal = QIn.Q / Q_nominal;
  annotation(Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10}), graphics = {Rectangle(visible = true, fillColor = {255, 0, 255}, fillPattern = FillPattern.Solid, extent = {{-50, -30}, {50, 30}})}));
end Terminal;
