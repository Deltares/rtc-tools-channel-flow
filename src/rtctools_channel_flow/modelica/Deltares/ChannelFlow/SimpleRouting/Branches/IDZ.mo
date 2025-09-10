within Deltares.ChannelFlow.SimpleRouting.Branches;

block IDZ
  import SI = Modelica.SIunits;
  extends Internal.PartialIDZ(H(min = H_b));

  // Bed level
  parameter SI.Position H_b;
  parameter SI.Time Delay_in_hour;

  annotation(Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10})));
end IDZ;
