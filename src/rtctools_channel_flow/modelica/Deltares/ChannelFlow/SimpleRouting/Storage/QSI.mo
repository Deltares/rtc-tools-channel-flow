within Deltares.ChannelFlow.SimpleRouting.Storage;

block QSI
  extends Deltares.ChannelFlow.Internal.QSI(QIn.Q(nominal=Q_nominal));
  extends Deltares.ChannelFlow.Internal.QForcing(QForcing(each nominal=Q_nominal));
  extends Deltares.ChannelFlow.Internal.Volume;
  parameter Modelica.Units.SI.VolumeFlowRate Q_nominal = 1.0;
equation
  // Mass balance
  der(V) / Q_nominal = (QIn.Q + sum(QForcing)) / Q_nominal;
end QSI;
