within Deltares.ChannelFlow.SimpleRouting.BoundaryConditions;

block Inflow
  import SI = Modelica.Units.SI;
  extends Deltares.ChannelFlow.Internal.QSO(QOut.Q(nominal=Q_nominal));
  // Inputs
  input SI.VolumeFlowRate Q(nominal=Q_nominal);

  parameter SI.VolumeFlowRate Q_nominal = 1.0;
equation
  QOut.Q / Q_nominal = Q / Q_nominal;
  annotation(Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10}), graphics = {Polygon(visible = true, fillColor = {255, 0, 255}, fillPattern = FillPattern.Solid, points = {{0, 50}, {-15, 15}, {-50, 0}, {-15, -15}, {0, -50}, {15, -15}, {50, 0}, {15, 15}})}));
end Inflow;
