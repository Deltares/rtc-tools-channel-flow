within Deltares.ChannelFlow.SimpleRouting.Storage; 

block Storage "DEPRECATED, use Branches.Integrator instead"
  import SI = Modelica.Units.SI;
  extends Deltares.ChannelFlow.Internal.QSISO(QIn.Q(nominal=Q_nominal), QOut.Q(nominal=Q_nominal));
  extends Deltares.ChannelFlow.Internal.QForcing(QForcing(each nominal = Q_nominal));
  // Inputs
  input SI.VolumeFlowRate Q_release(nominal=Q_nominal);
  // States
  SI.Volume V(min=0, nominal = 1e6);

  parameter SI.VolumeFlowRate Q_nominal = 1.0;
equation
  // Mass balance
  der(V) / Q_nominal = (QIn.Q - QOut.Q + sum(QForcing)) / Q_nominal;
  // Outflow equals release
  QOut.Q / Q_nominal = Q_release / Q_nominal;
  annotation(Icon(coordinateSystem(initialScale = 0.1, grid = {10, 10}), graphics = {Rectangle(fillColor = {255, 0, 0}, fillPattern = FillPattern.Solid, extent = {{-50, 50}, {50, -50}})}));
end Storage;
