within Deltares.ChannelFlow.Hydraulic.Branches;

model IDZ
  extends Deltares.ChannelFlow.Hydraulic.Branches.Internal.PartialIDZ;

  // Bed level
  parameter Modelica.Units.SI.Position H_b_up;
  parameter Modelica.Units.SI.Position H_b_down;
  parameter Modelica.Units.SI.Position length;
  parameter Modelica.Units.SI.Position width;
  parameter Modelica.Units.SI.VolumeFlowRate Q_nominal;
  parameter Modelica.Units.SI.Position H_nominal;
  parameter Modelica.Units.SI.Position friction_coefficient;
  parameter Modelica.Units.SI.Position side_slope;
  parameter Integer n_level_nodes = 2;

  annotation(Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10})));
end IDZ;
