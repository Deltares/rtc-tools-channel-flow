within Deltares.ChannelFlow.SimpleRouting.Storage;

block QSI
  extends Deltares.ChannelFlow.Internal.QSI;
  extends Deltares.ChannelFlow.Internal.QForcing;
  extends Deltares.ChannelFlow.Internal.Volume;
equation
  // Mass balance
  der(V) = QIn.Q + sum(QForcing);
end QSI;
