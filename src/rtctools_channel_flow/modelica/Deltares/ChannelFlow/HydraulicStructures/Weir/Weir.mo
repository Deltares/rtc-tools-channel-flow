within Deltares.ChannelFlow.HydraulicStructures.Weir;

model Weir
  import SI = Modelica.Units.SI;
  extends Deltares.ChannelFlow.Internal.HQTwoPort;

  // Inputs
  input SI.VolumeFlowRate Q;

  // Parameters
  parameter SI.Length width "Width of the weir";
  parameter SI.VolumeFlowRate q_min "Minimum flow of the weir; has to be positive";
  parameter SI.VolumeFlowRate q_max "Maximum flow of the weir. Should be as low as possible.";
  parameter SI.Length hw_min "Minimum height of the weir";
  parameter SI.Length hw_max "Maximum height of the weir";
  parameter Real weir_coef=0.61 "Weir discharge coefficient";
equation
  HQUp.Q + HQDown.Q = 0; // Negative comes in, positive out, so in a branch positive goes in
  HQUp.Q = Q;
  annotation(Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10}), graphics = {Polygon(visible = true, origin = {0, -16.667}, fillColor = {255, 128, 0}, fillPattern = FillPattern.Solid, lineThickness = 2, points = {{0, 66.667}, {-50, -33.333}, {50, -33.333}})}), Diagram(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10})));
end Weir;
