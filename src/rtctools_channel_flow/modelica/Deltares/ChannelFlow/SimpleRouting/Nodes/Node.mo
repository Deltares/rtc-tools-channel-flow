within Deltares.ChannelFlow.SimpleRouting.Nodes;

block Node "Block with multiple inflows and multiple outflows and forcing, where allocation is based on explicitly specified outflows."
  import SI = Modelica.Units.SI;
  extends Internal.PartialNode(redeclare parameter Integer nout(min = 1) = 1);
  extends Deltares.ChannelFlow.Internal.QForcing(QForcing(each nominal=Q_nominal));
  input SI.VolumeFlowRate QOut_control[nout - 1](each nominal=Q_nominal);
equation
  QInSum / Q_nominal = sum(QIn.Q) / Q_nominal;
  QOutSum / Q_nominal = (QInSum + sum(QForcing)) / Q_nominal;
  for i in 1:nout - 1 loop
    QOut[i].Q / Q_nominal = QOut_control[i] / Q_nominal;
  end for;
  QOut[nout].Q / Q_nominal = (QOutSum - sum(QOut_control[1:nout - 1])) / Q_nominal;
  annotation(Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10}), graphics = {Polygon(visible = true, fillColor = {255, 170, 0}, fillPattern = FillPattern.Solid, points = {{0, 50}, {-30, 40}, {30, -40}, {0, -50}, {-30, -40}, {30, 40}}), Polygon(visible = true, fillColor = {255, 0, 0}, fillPattern = FillPattern.Solid, points = {{-50, 0}, {-40, 30}, {-30, 40}, {30, -40}, {40, -30}, {50, 0}, {40, 30}, {30, 40}, {-30, -40}, {-40, -30}})}), Diagram(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10})));
end Node;
