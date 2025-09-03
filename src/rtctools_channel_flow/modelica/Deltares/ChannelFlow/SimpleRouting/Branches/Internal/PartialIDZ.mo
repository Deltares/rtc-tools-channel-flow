within Deltares.ChannelFlow.SimpleRouting.Branches.Internal;

partial block PartialIDZ
  extends Deltares.ChannelFlow.Internal.HQTwoPort;
  extends Deltares.ChannelFlow.Internal.QForcing;
  extends Deltares.ChannelFlow.Internal.QLateral;
  //extends Deltares.ChannelFlow.Internal.Reservoir;
  input Real Q_upstream_delayed;
  input Real Q_downstream_delayed;

  // States
  Modelica.SIunits.Position H;
  parameter Real p21;
  parameter Real p22;
  parameter Real p11;
  parameter Real p12;
  parameter SI.Area Au;
  parameter SI.Area Ad;

equation
  // Water level
  H = HQUp.H;

  der(HQDown.H) = Q_upstream_delayed / Ad + sum(QForcing) / Ad + sum(QLateral.Q) / Ad + p21 * der(Q_upstream_delayed)+ HQDown.Q / Ad + p22 * der(HQDown.Q);
  der(HQUp.H) =   HQUp.Q / Au + sum(QForcing) / Au + sum(QLateral.Q) / Au + p11 * der(HQUp.Q)+ Q_downstream_delayed / Au + p12 * der(Q_downstream_delayed);

end PartialIDZ;