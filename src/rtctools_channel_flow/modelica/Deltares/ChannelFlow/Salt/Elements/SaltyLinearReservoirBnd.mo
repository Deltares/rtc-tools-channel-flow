within Deltares.ChannelFlow.Salt.Elements;

model SaltyLinearReservoirBnd
  /*
  This block is designed to be used together with the "salt_simulation_mixin" to calculate dispersive and advective transport
   between salty reservoir elements, do not use in optimization.
  */
  import SI = Modelica.Units.SI;
  extends SaltyPartialReservoirBnd(H(min = H_b));
  
  parameter SI.Area A;
  parameter SI.Position H_b; // Bed level

equation
  
  V = A * (H - H_b);
  annotation(Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10})));
end SaltyLinearReservoirBnd;
