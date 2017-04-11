within Deltares.ChannelFlow.Hydraulic.Structures;

model Pump "Pump"
  extends Deltares.ChannelFlow.Internal.HQCMTwoPort;
  input Modelica.SIunits.VolumeFlowRate Q;
  parameter Modelica.SIunits.MassFlowRate Q_nominal = 1;
  parameter Modelica.SIunits.Density[HQCMUp.NOS] C_nominal = fill(1,HQCMUp.NOS);
  parameter Real theta  = 1;
equation
  HQCMUp.Q + HQCMDown.Q = 0;
  HQCMUp.Q = Q;
  
  HQCMUp.M = -HQCMDown.M;
  //Z depends on which direction the flow is, this decouples the concentration on both sides of the pump.
  //z=Q*C, this equation is linearized.
  if(Q > 0) then
    HQCMUp.M = theta * HQCMUp.C * Q + (1 - theta) * (Q_nominal * C_nominal + C_nominal * (Q - Q_nominal) + Q_nominal * (HQCMUp.C - C_nominal));
  else 
    HQCMUp.M = theta * HQCMDown.C * Q + (1 - theta) * (Q_nominal * C_nominal + C_nominal * (Q - Q_nominal) + Q_nominal * (HQCMDown.C - C_nominal));
  end if; 
  annotation(Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10}), graphics = {Polygon(visible = true, origin = {0, -16.667}, fillColor = {255, 128, 0}, fillPattern = FillPattern.Solid, lineThickness = 2, points = {{0, 66.667}, {-50, -33.333}, {50, -33.333}})}), Diagram(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10})));
end Pump;
