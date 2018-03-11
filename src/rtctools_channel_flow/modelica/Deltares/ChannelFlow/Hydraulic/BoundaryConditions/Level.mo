within Deltares.ChannelFlow.Hydraulic.BoundaryConditions;

model Level "Defines absolute water level"
  extends Deltares.ChannelFlow.Internal.HQOnePort;
  input Modelica.SIunits.Position H;
  input Modelica.SIunits.Density C[HQ.medium.n_substances];
equation
  HQ.H = H;
  HQ.C = C;
  annotation(__Wolfram(itemFlippingEnabled = true), Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10}), graphics = {Rectangle(visible = true, fillColor = {255, 0, 255}, fillPattern = FillPattern.Solid, extent = {{-50, -50}, {50, 50}})}));
end Level;
