within Deltares.ChannelFlow.Hydraulic.Reservoir.Internal;

partial model PartialReservoir
  import SI = Modelica.SIunits;
  extends Deltares.ChannelFlow.Internal.HQTwoPort;
  extends Deltares.ChannelFlow.Internal.QForcing;
  // Inputs
  input SI.VolumeFlowRate Q_turbine;
  input SI.VolumeFlowRate Q_spill;
  // States
  SI.Position H;
  SI.Volume V(min=0, nominal = 1e6);
equation
  // Water level
  H = HQUp.H;
  // Mass balance
  der(V) = HQUp.Q + HQDown.Q + sum(QForcing);
  // Split outflow between turbine and spill flow
  HQDown.Q + Q_turbine + Q_spill = 0.0;
  annotation(Icon(coordinateSystem( initialScale = 0.1, grid = {10, 10}), graphics = {Polygon( fillColor = {0, 255, 255}, fillPattern = FillPattern.Solid, points = {{40, 50}, {-45, 0}, {40, -50}, {40, 50}})}));
end PartialReservoir;