within Deltares.ChannelFlow.Hydraulic.Branches.Internal;

partial model PartialIDZ
  extends Deltares.ChannelFlow.Internal.HQTwoPort;
  extends Deltares.ChannelFlow.Internal.QForcing;
  extends Deltares.ChannelFlow.Internal.QLateral;
  extends Deltares.ChannelFlow.Internal.Reservoir;
  input Real Q_upstream_delayed;
  input Real Q_downstream_delayed;
  //output Real q_up_der;

  // States
  SI.Position[2] H(min = H_b_down);
  SI.VolumeFlowRate[2] Q; 
  parameter Real p21;
  parameter Real p22;
  parameter Real p11;
  parameter Real p12;
  parameter SI.Area Au;
  parameter SI.Area Ad;

equation
  // Water level
  H[1] = HQUp.H;
  H[2] = HQDown.H;
  Q[1] = HQUp.Q;
  Q[2] = HQDown.Q;
  
  // Mass balance
  //der(V) = HQUp.Q + HQDown.Q + sum(QForcing) + sum(QLateral.Q);
  //idz_zero = 0.001;
  //der(HQDown.Q) = q_up_der;
  //der(HQDown.Q) =   q_down_der;

  der(HQDown.H) = Q_upstream_delayed / Ad + sum(QForcing) / Ad + sum(QLateral.Q) / Ad + p21 * der(Q_upstream_delayed)+ HQDown.Q / Ad + p22 * der(HQDown.Q);
  der(HQUp.H) =   HQUp.Q / Au + sum(QForcing) / Au + sum(QLateral.Q) / Au + p11 * der(HQUp.Q)+ Q_downstream_delayed / Au + p12 * der(Q_downstream_delayed);
  Q_upstream_delayed = HQUp.Q;
  Q_downstream_delayed = HQDown.Q;

  //Q_upstream_delayed + HQDown.Q = 0;

  // Split outflow between turbine and spill flow
  HQDown.Q + Q_turbine + Q_spill = 0.0;
end PartialIDZ;