within Deltares.ChannelFlow.SimpleRouting.Reservoir;

block Reservoir
  import SI = Modelica.SIunits;
  extends Deltares.ChannelFlow.Internal.QSISO;
  extends Deltares.ChannelFlow.Internal.QForcing;
  // Inputs
  input SI.VolumeFlowRate Q_turbine;
  input SI.VolumeFlowRate Q_spill;
  // States
  SI.Volume V(nominal = 1e6);
equation
  // Mass balance
  der(V) = QIn.Q - QOut.Q + sum(QForcing);
  // Split outflow between turbine and spill flow
  QOut.Q = Q_turbine + Q_spill;
  annotation(Icon(coordinateSystem( initialScale = 0.1, grid = {10, 10}), graphics = {Polygon(fillColor = {0, 255, 255}, fillPattern = FillPattern.Solid, points = {{40, 50}, {-45, 0}, {40, -50}, {40, 50}, {40, 50}})}));
end Reservoir;