within Deltares.ChannelFlow.Salt.Elements;

model NodeSalty 
  /*
  This block is designed to be used together with the "salt_simulation_mixin" to calculate dispersive and advective transport
   between salty reservoir elements, do not use in optimization.
  */
  extends Deltares.ChannelFlow.Internal.HQTwoPort;
  extends QMForcing;

equation
  
  HQUp.Q + HQDown.Q + sum(QForcing) = 0;
  HQUp.M[1] + HQDown.M[1] + sum(MForcing)  = 0; //Mass balance, can be used for one substance, like salt.
    
  HQDown.H = HQUp.H;
  HQUp.C = HQDown.C;

  annotation(Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10}), graphics = {Polygon(visible = true, fillColor = {255, 170, 0}, fillPattern = FillPattern.Solid, points = {{0, 50}, {-30, 40}, {30, -40}, {0, -50}, {-30, -40}, {30, 40}}), Polygon(visible = true, fillColor = {255, 0, 0}, fillPattern = FillPattern.Solid, points = {{-50, 0}, {-40, 30}, {-30, 40}, {30, -40}, {40, -30}, {50, 0}, {40, 30}, {30, 40}, {-30, -40}, {-40, -30}})}));
end NodeSalty;
