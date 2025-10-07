within Deltares.ChannelFlow.Hydraulic.BoundaryConditions;

model Level "Defines absolute water level"
  /*
  Note: The default medium is FreshWater.
  To use a different medium, decalre the choice in your model file, for example
  replaceable package MyMedium = Deltares.ChannelFlow.Media.SalineWater;
  Pass this as an argument to the Level block (redeclare package medium=MyMedium)
  */
  extends Deltares.ChannelFlow.Internal.HQOnePort;
  input Modelica.Units.SI.Position H;
  input Modelica.Units.SI.Density C[medium.n_substances];
equation
  HQ.H = H;
  HQ.C = C;
  annotation(__Wolfram(itemFlippingEnabled = true), Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10}), graphics = {Rectangle(visible = true, fillColor = {255, 0, 255}, fillPattern = FillPattern.Solid, extent = {{-50, -50}, {50, 50}})}));
end Level;
