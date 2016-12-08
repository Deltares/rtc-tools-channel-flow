within Deltares.ChannelFlow.SimpleRouting.Storage;

block QSO
  import SI = Modelica.SIunits;
  extends Deltares.ChannelFlow.Internal.QSO;
  extends Deltares.ChannelFlow.Internal.QForcing;
  // Inputs
  input SI.VolumeFlowRate QOut_control;
  // States
  SI.Volume V(min=0, nominal = 1e6);
equation
  // Mass balance
  der(V) = - QOut.Q + sum(QForcing);
  // Outflow equals release
  QOut.Q = QOut_control;
  annotation(Icon(coordinateSystem(initialScale = 0.1, grid = {10, 10}), graphics = {Rectangle(fillColor = {255, 0, 0}, fillPattern = FillPattern.Solid, extent = {{-50, 50}, {50, -50}})}));
end QSO;
