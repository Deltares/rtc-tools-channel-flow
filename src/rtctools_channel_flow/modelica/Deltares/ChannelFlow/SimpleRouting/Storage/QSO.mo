within Deltares.ChannelFlow.SimpleRouting.Storage;

block QSO
  import Modelica.Units.SI;
  extends Deltares.ChannelFlow.Internal.QSO;
  extends Deltares.ChannelFlow.Internal.QForcing;
  extends Deltares.ChannelFlow.Internal.Volume;
  // Inputs
  input SI.VolumeFlowRate QOut_control;
equation
// Mass balance
  der(V) = -QOut.Q + sum(QForcing);
// Outflow equals release
  QOut.Q = QOut_control;
end QSO;