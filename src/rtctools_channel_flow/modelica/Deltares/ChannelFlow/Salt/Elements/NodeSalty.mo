within Deltares.ChannelFlow.Salt.Elements;

block NodeSalty "Block with multiple inflows and multiple outflows, where allocation is based on explicitly specified outflows, including a port for a reservoir, includes concentration"

  extends Deltares.ChannelFlow.Internal.HQTwoPort;
  extends QMForcing;

  // States
  //SI.Concentration C(nominal = 0.001);
  //SI.Height H(nominal = 0.001);



equation
  //HQUp.Q + HQDown.Q  = 0;
  //HQUp.M + HQDown.M   = 0;
  
  HQUp.Q + HQDown.Q + sum(QForcing) = 0;
  HQUp.M + HQDown.M + sum(MForcing)  = 0;
  
  
  HQDown.H = HQUp.H;
  HQUp.C = HQDown.C;
  //HQUp.Q =0;
  //Q = HQDown.Q;
  //HQUp.M =0;
  //C = 0;
  //H = 0;

  annotation(Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10}), graphics = {Polygon(visible = true, fillColor = {255, 170, 0}, fillPattern = FillPattern.Solid, points = {{0, 50}, {-30, 40}, {30, -40}, {0, -50}, {-30, -40}, {30, 40}}), Polygon(visible = true, fillColor = {255, 0, 0}, fillPattern = FillPattern.Solid, points = {{-50, 0}, {-40, 30}, {-30, 40}, {30, -40}, {40, -30}, {50, 0}, {40, 30}, {30, 40}, {-30, -40}, {-40, -30}})}));
end NodeSalty;
