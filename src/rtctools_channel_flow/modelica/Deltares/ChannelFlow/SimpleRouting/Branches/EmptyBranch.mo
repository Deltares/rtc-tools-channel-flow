within Deltares.ChannelFlow.SimpleRouting.Branches;

model EmptyBranch "Branch whose delays can be modified within the python src"
  import SI = Modelica.Units.SI;
  extends Deltares.ChannelFlow.Internal.QSISO(QIn.Q(nominal=Q_nominal), QOut.Q(nominal=Q_nominal));
  // Nominal values for scaling
  parameter SI.VolumeFlowRate Q_nominal = 1.0;
equation
  annotation(Icon(coordinateSystem( initialScale = 0.1, grid = {10, 10}), graphics = {Line(points = {{-50, 0}, {50, 0}})}));
  end EmptyBranch;
