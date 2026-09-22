within Deltares.ChannelFlow.Hydraulic.Branches;

model ID
  import SI = Modelica.Units.SI;
  extends Deltares.ChannelFlow.Internal.HQTwoPort;
  extends Deltares.ChannelFlow.Internal.QForcing;
  extends Deltares.ChannelFlow.Internal.QLateral;
  extends Deltares.ChannelFlow.Internal.Reservoir;

  input Real Q_upstream_delayed;
  input Real Q_downstream_delayed;
  parameter SI.Duration Delay_in_hour;

  // States
  SI.Position[2] H;
  SI.VolumeFlowRate[2] Q;
  parameter SI.Area Ad;

equation
  // Water level
  H[1] = HQUp.H;
  H[2] = HQDown.H;
  Q[1] = HQUp.Q;
  Q[2] = HQDown.Q;

  der(HQDown.H) = Q_upstream_delayed / Ad + sum(QForcing) / Ad + sum(QLateral.Q) / Ad + HQDown.Q / Ad ;
  der(HQUp.H) =   HQUp.Q / Ad + sum(QForcing) / Ad + sum(QLateral.Q) / Ad + Q_downstream_delayed / Ad ;

  Q_downstream_delayed = delay(HQDown.Q, Delay_in_hour);
  Q_upstream_delayed = delay(HQUp.Q, Delay_in_hour);


  annotation(Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10})), Documentation(info="
  <html>
  <p>
  This block represents an Integrator Delay (ID) reach, modelling the propagation
  and storage of water between an upstream and downstream location.
  </p>
  <p>
  Water levels are represented at both ends of the reach. The rate of change of
  each water level is determined by the local discharge, external forcings, lateral
  inflows, and a delayed discharge originating from the opposite boundary.
  The delay represents the travel time of a flow wave through the reach.
  </p>
  <p>
  Delayed discharges are computed using the Modelica <code>delay()</code> operator:
  </p>
  <math>
  Q_{up}^{delayed} = delay(Q_{up}, T_d)
  </math>
  <math>
  Q_{down}^{delayed} = delay(Q_{down}, T_d)
  </math>
  <p>
  where <code>T_d</code> is the travel time through the reach.
  </p>
  <p>
  The resulting delayed flows act as boundary conditions for the opposite side
  of the reach, allowing the model to represent transport and attenuation effects
  while maintaining a simple storage-based formulation.
  </p>
  </html>"));
end ID;
