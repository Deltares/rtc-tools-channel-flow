within Deltares.ChannelFlow.Hydraulic.Reservoir.Internal;

partial model PartialHomotopicVolume
  import SI = Modelica.Units.SI;
  extends Deltares.ChannelFlow.Hydraulic.Reservoir.Internal.PartialReservoir;

  parameter SI.Area A;
  // Bed level
  parameter SI.Position H_b;
  // Homotopy parameter
  parameter Real theta;
  // Water level polynomial coefficients, to be specified as parameters in the model
  parameter Real Hc0 = 0.0;
  parameter Real Hc1 = 0.0;
  parameter Real Hc2 = 0.0;
  parameter Real Hc3 = 0.0;
  parameter Real Hc4 = 0.0;
equation
  // Volume - forebay relation
  V / A = ((1 - theta) * A * (H - H_b) + theta * (Hc0 + Hc1*H + Hc2*H^2 + Hc3*H^3 + Hc4*H^4)) / A;
end PartialHomotopicVolume;
