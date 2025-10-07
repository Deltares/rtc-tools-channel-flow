within Deltares.ChannelFlow.SimpleRouting.Branches;

block Integrator
  import SI = Modelica.Units.SI;
  extends Deltares.ChannelFlow.Internal.QSISO(QIn.Q(nominal=Q_nominal), QOut.Q(nominal=Q_nominal));
  extends Deltares.ChannelFlow.Internal.QForcing(QForcing(each nominal = Q_nominal));
  extends Deltares.ChannelFlow.Internal.QLateral(QLateral.Q(each nominal=Q_nominal));
  extends Deltares.ChannelFlow.Internal.Volume;
  // Inputs
  input SI.VolumeFlowRate QOut_control(nominal=Q_nominal);
  parameter SI.VolumeFlowRate Q_nominal = 1.0;
equation
  // Mass balance
  der(V) / Q_nominal = (QIn.Q - QOut.Q + sum(QForcing) + sum(QLateral.Q)) / Q_nominal;

  // Outflow equals release
  QOut.Q / Q_nominal = QOut_control / Q_nominal;
end Integrator;
