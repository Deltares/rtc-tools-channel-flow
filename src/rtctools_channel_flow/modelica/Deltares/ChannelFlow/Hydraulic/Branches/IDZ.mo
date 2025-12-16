within Deltares.ChannelFlow.Hydraulic.Branches;

model IDZ
  import SI = Modelica.Units.SI;
  extends PartialIDZ;

  // Bed level
  parameter SI.Position H_b_up;
  parameter SI.Position H_b_down;
  parameter SI.Time Delay_in_hour = 0.0;
  parameter SI.Position length;
  parameter SI.Position width;
  parameter SI.VolumeFlowRate Q_nominal;
  parameter SI.Position H_nominal;
  parameter SI.Position friction_coefficient;
  parameter SI.Position side_slope;
  //parameter n_level_nodes = 2;

  parameter Integer n_level_nodes = 2;

  annotation(Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10})));
end IDZ;
