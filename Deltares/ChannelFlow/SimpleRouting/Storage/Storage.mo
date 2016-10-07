within Deltares.ChannelFlow.SimpleRouting.Storage;

block Storage
  import SI = Modelica.SIunits;
  extends Deltares.ChannelFlow.Internal.QSISO;
  extends Deltares.ChannelFlow.Internal.QForcing;
  // Inputs
  input SI.VolumeFlowRate Q_release;
  // States
  SI.Volume V(nominal = 1e6);
equation
  // Mass balance
  der(V) = QIn.Q - QOut.Q + sum(QForcing);
  // Outflow equals release
  QOut.Q = Q_release;
  annotation(Icon(coordinateSystem(initialScale = 0.1, grid = {10, 10}), graphics = {Rectangle(fillColor = {255, 0, 0}, fillPattern = FillPattern.Solid, extent = {{-50, 50}, {50, -50}})}));
end Storage;
