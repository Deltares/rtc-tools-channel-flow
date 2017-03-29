within Deltares.ChannelFlow.Hydraulic.Structures;

model PumpZ
  extends Deltares.ChannelFlow.Hydraulic.Structures.Pump(redeclare connector HQPort = Deltares.ChannelFlow.Interfaces.HQCMPort);
  parameter Modelica.SIunits.MassFlowRate Q_nominal = 1;
  parameter Modelica.SIunits.Density C_nominal = 1;
  parameter Real theta;
equation
  HQUp.M + HQDown.M = 0;
  //Z depends on which direction the flow is, this decouples the concentration on both sides of the pump.
  //z=Q*C, this equation is linearized.
  if(Q > 0) then
    HQUp.M = theta * HQUp.C * Q + (1 - theta) * (Q_nominal * C_nominal + C_nominal * (Q - Q_nominal) + Q_nominal * (HQUp.C - C_nominal));
  else
    HQUp.M = theta * HQDown.C * Q + (1 - theta) * (-Q_nominal * C_nominal + C_nominal * (Q + Q_nominal) - Q_nominal * (HQDown.C - C_nominal));
  end if;                               
end PumpZ;
