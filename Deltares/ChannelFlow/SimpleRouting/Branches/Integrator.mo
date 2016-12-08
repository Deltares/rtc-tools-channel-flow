within Deltares.ChannelFlow.SimpleRouting.Branches;

block Integrator
  import SI = Modelica.SIunits;
  extends Deltares.ChannelFlow.Internal.QSISO;
  extends Deltares.ChannelFlow.Internal.QForcing;
  extends Deltares.ChannelFlow.Internal.QLateral;
  // Inputs
  input SI.VolumeFlowRate QOut_control;
  // States
  SI.Volume V(min=0, nominal = 1e6);
equation
  // Mass balance
  der(V) = QIn.Q - QOut.Q + sum(QForcing) + sum(QLateral.Q);
  // Outflow equals release
  QOut.Q = QOut_control;
  annotation(Icon(coordinateSystem(initialScale = 0.1, grid = {10, 10}), graphics = {Rectangle(fillColor = {255, 0, 0}, fillPattern = FillPattern.Solid, extent = {{-50, 50}, {50, -50}})}));
end Integrator;
