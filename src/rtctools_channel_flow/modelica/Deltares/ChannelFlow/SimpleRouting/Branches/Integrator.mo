within Deltares.ChannelFlow.SimpleRouting.Branches;

block Integrator
  import SI = Modelica.SIunits;
  extends Deltares.ChannelFlow.Internal.QSISO;
  extends Deltares.ChannelFlow.Internal.QForcing;
  extends Deltares.ChannelFlow.Internal.QLateral;
  extends Deltares.ChannelFlow.Internal.Volume;
  // Inputs
  input SI.VolumeFlowRate QOut_control;
equation
  // Mass balance
  der(V) = QIn.Q - QOut.Q + sum(QForcing) + sum(QLateral.Q);
  // Outflow equals release
  QOut.Q = QOut_control;
end Integrator;
