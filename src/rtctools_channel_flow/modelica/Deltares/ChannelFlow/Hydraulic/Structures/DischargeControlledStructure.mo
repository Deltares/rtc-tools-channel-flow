within Deltares.ChannelFlow.Hydraulic.Structures;

model DischargeControlledStructure "DischargeControlledStructure"
  extends Deltares.ChannelFlow.Internal.HQTwoPort;
  function smooth_switch = Deltares.ChannelFlow.Internal.Functions.SmoothSwitch;
  input Modelica.SIunits.VolumeFlowRate Q;
  // Homotopy parameter
  parameter Real theta = 1.0;
  // Nominal values used in linearization
  parameter Modelica.SIunits.MassFlowRate Q_nominal = 1;
  parameter Modelica.SIunits.Density C_nominal[HQUp.medium.n_substances] = fill(1e-3, HQUp.medium.n_substances);
equation
  // Water
  HQUp.Q + HQDown.Q = 0;
  HQUp.Q = Q;
  // Substances
  HQUp.M = -HQDown.M;
  // Z depends on which direction the flow is, this decouples the concentration on both sides of the pump.
  // Z=Q*C, this equation is linearized.
  HQUp.M = theta * (smooth_switch(Q) * HQUp.C * Q + (1 - smooth_switch(Q)) * HQDown.C * Q) + (1 - theta) * (Q_nominal * C_nominal + C_nominal * (Q - Q_nominal) + Q_nominal * ((if Q_nominal > 0 then HQUp.C else HQDown.C) - C_nominal));
  annotation(Icon(coordinateSystem( initialScale = 0.1, grid = {10, 10}), graphics = {Polygon(origin = {0, -16.67}, fillColor = {255, 128, 0}, fillPattern = FillPattern.Solid, points = {{0, 66.667}, {-50, -33.333}, {50, -33.333}, {0, 66.667}})}), Diagram(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10})));
end DischargeControlledStructure;