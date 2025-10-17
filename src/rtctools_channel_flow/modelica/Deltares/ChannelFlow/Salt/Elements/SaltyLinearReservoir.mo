within Deltares.ChannelFlow.Salt.Elements;

model SaltyLinearReservoir
  /*
  This block is designed to be used together with the "salt_simulation_mixin" to calculate dispersive and advective transport
   between salty reservoir elements, do not user in optimization.
  */
  import SI = Modelica.Units.SI;
  extends SaltyPartialReservoir(H(min = H_b));
  // Parameters
  parameter SI.Area A;
  // Bed level
  parameter SI.Position H_b;
equation
  // Volume - forebay relation
  V = A * (H - H_b);
  annotation(Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10})));
end SaltyLinearReservoir;
