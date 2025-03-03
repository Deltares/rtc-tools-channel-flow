within Deltares.ChannelFlow.Hydraulic.Reservoir.Internal;

partial model PartialHomotopicVolume
// this partial model adds a nonlinear non convex relation between volume and water level to the reservoir. 
// the relation between water level and volume is specified as a fourth order polynomial. 
// the optimization starts with a linearized simplification, here the volume is divided by a area A to obtain the water level. 
// during the homotopy optimization this linear equation is gradually bended towards the fourth order polynomial with the help of theta. 
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
