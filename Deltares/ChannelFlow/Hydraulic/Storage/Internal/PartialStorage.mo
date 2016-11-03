within Deltares.ChannelFlow.Hydraulic.Storage.Internal;

partial model PartialStorage
  extends Deltares.ChannelFlow.Internal.HQOnePort;
  extends Deltares.ChannelFlow.Internal.QForcing;
  Modelica.SIunits.Volume V(min=0, nominal = 1e6);
equation
  der(V) = HQ.Q + sum(QForcing);
  annotation(Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10}), graphics = {Rectangle(visible = true, fillColor = {255, 0, 0}, fillPattern = FillPattern.Solid, extent = {{-50, -50}, {50, 50}})}), Diagram(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10})));
end PartialStorage;
