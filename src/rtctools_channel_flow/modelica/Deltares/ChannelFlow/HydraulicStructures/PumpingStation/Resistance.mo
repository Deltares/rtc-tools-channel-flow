within Deltares.ChannelFlow.HydraulicStructures.PumpingStation;

// TODO: Negative flows (from down to up) are not supported. Do we want to support them?
model Resistance "Quadratic resistance of form dH=C*Q^2"
  extends Deltares.ChannelFlow.Internal.HQTwoPort;

  parameter Real C = 0.0;

  // Head loss
  input Modelica.Units.SI.Distance dH;
equation
  // Head
  HQDown.H = HQUp.H - dH;

  // Discharge
  HQUp.Q + HQDown.Q = 0;

  // Substances
  HQUp.M = -HQDown.M;

  // TODO: Annotation / pretty picture. Currently inheriting TwoPort.
end Resistance;
