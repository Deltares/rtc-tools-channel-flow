within Deltares.ChannelFlow.SimpleRouting.Storage;

block QSI
  import SI = Modelica.SIunits;
  extends Deltares.ChannelFlow.Internal.QSI;
  extends Deltares.ChannelFlow.Internal.QForcing;
  // States
  SI.Volume V(min=0, nominal = 1e6);
equation
  // Mass balance
  der(V) = QIn.Q + sum(QForcing);
  annotation(Icon(coordinateSystem(initialScale = 0.1, grid = {10, 10}), graphics = {Rectangle(fillColor = {255, 0, 0}, fillPattern = FillPattern.Solid, extent = {{-50, 50}, {50, -50}})}));
end QSI;
