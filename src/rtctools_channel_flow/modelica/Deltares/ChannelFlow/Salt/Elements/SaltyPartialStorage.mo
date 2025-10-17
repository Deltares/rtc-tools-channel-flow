within Deltares.ChannelFlow.Salt.Elements;

partial model SaltyPartialStorage
  import SI = Modelica.Units.SI;
  extends Deltares.ChannelFlow.Internal.HQOnePort(HQ.Q(each nominal = Q_nominal), HQ.M(each nominal =Q_nominal - C_nominal));
  extends Deltares.ChannelFlow.Internal.QForcing(QForcing(each nominal = Q_nominal));
  extends Deltares.ChannelFlow.Internal.Volume;
  // Homotopy parameter
  parameter Real theta = 1.0;
  // Nominal values used in linearization
  parameter SI.Volume V_nominal;
  parameter SI.Density C_nominal = 1e-3;
  // Nominal values for scaling
  parameter SI.VolumeFlowRate Q_nominal = 1.0;
equation
  der(V) / Q_nominal = (HQ.Q + sum(QForcing)) / Q_nominal;
  HQ.M / (Q_nominal * C_nominal) = (theta * der(V * HQ.C) + (1 - theta) * Q_nominal * der(HQ.C))  / (Q_nominal * C_nominal);
end SaltyPartialStorage;
