within Deltares.ChannelFlow.Hydraulic.Reservoir;

model Linear
  import SI = Modelica.SIunits;
  extends Internal.PartialReservoir(H(min = H_b));
  // Parameters
  parameter SI.Area A;
  // Bed level
  parameter SI.Position H_b;
equation
  // Volume - forebay relation
  V = A * (H - H_b);
  annotation(Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10})));
end Linear;
