within Deltares.ChannelFlow.Hydraulic.Storage.Internal;

partial model PartialStorage
  extends Deltares.ChannelFlow.Internal.HQOnePort;
  extends Deltares.ChannelFlow.Internal.QForcing;
  extends Deltares.ChannelFlow.Internal.Volume;
equation
  der(V) = HQ.Q + sum(QForcing);
end PartialStorage;
