within Deltares.ChannelFlow.Hydraulic.Reservoir.Internal;

partial model PartialHomotopicVolume
  import SI = Modelica.Units.SI;
  extends Deltares.ChannelFlow.Internal.HQTwoPort(HQDown.Q(nominal = Q_nominal));
  extends Deltares.ChannelFlow.Internal.QForcing;
  extends Deltares.ChannelFlow.Internal.QLateral;
  extends Deltares.ChannelFlow.Internal.Reservoir;
  // Parameters
  parameter Real Q_nominal=1.0;
  parameter SI.Area A;
  // Water level
  Modelica.Units.SI.Position H(min = H_b) = HQUp.H;
  // Bed level
  parameter SI.Position H_b;
  // Homotopy parameter
  parameter Real theta;
  // Water level polynomial coefficients
  parameter Real Hc0 = 0.0;
  parameter Real Hc1 = 0.0;
  parameter Real Hc2 = 0.0;
  parameter Real Hc3 = 0.0;
  parameter Real Hc4 = 0.0;
equation
  // Mass balance
  der(V) / Q_nominal = (HQUp.Q + HQDown.Q + sum(QForcing) + sum(QLateral.Q)) / Q_nominal;
  // Volume - forebay relation
  V / A = ((1 - theta) * A * (H - H_b) + theta * (Hc0 + Hc1*H + Hc2*H^2 + Hc3*H^3 + Hc4*H^4)) / A;
  // Split outflow between turbine and spill flow
  (HQDown.Q + Q_turbine + Q_spill) / Q_nominal = 0.0;
end PartialHomotopicVolume;
