within Deltares.ChannelFlow.SimpleRouting.Branches;

model EmptyBranch "Branch whose delays can be modified within the python src"
  extends Deltares.ChannelFlow.Internal.QSISO;
equation
  annotation(Icon(coordinateSystem( initialScale = 0.1, grid = {10, 10}), graphics = {Line(points = {{-50, 0}, {50, 0}})}));
  end EmptyBranch;
