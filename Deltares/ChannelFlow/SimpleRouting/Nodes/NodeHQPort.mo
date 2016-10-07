within Deltares.ChannelFlow.SimpleRouting.Nodes;

block NodeHQPort "Block with multiple inflows and multiple outflows, where allocation is based on explicitly specified outflows, including a port for a reservoir"
  import SI = Modelica.SIunits;
  extends Internal.PartialNode;
  extends Deltares.ChannelFlow.Internal.HQOnePort;
  input SI.VolumeFlowRate QOut_control[nout];
  output SI.Position H;
equation
  QInSum = sum(QIn.Q);
  QOutSum = sum(QOut_control);
  for i in 1:nout loop
    QOut[i].Q = QOut_control[i];
  end for;
  HQ.Q + QInSum - QOutSum = 0.;
  HQ.H = H;
  annotation(Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10}), graphics = {Polygon(visible = true, fillColor = {255, 170, 0}, fillPattern = FillPattern.Solid, points = {{0, 50}, {-30, 40}, {30, -40}, {0, -50}, {-30, -40}, {30, 40}}), Polygon(visible = true, fillColor = {255, 0, 0}, fillPattern = FillPattern.Solid, points = {{-50, 0}, {-40, 30}, {-30, 40}, {30, -40}, {40, -30}, {50, 0}, {40, 30}, {30, 40}, {-30, -40}, {-40, -30}})}));
end NodeHQPort;
