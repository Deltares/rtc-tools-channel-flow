within Deltares.ChannelFlow.Hydraulic.BoundaryConditions;

model DischargeMpos
  extends Deltares.ChannelFlow.Hydraulic.BoundaryConditions.Discharge(redeclare connector HQPort = Deltares.ChannelFlow.Interfaces.HQCMPort);
  input Modelica.SIunits.Density[HQ.NOS] C;
equation
  if(Q>0) then
    HQ.M = -Q * C;
  else
    HQ.M = -Q * HQ.C;
  end if;
end DischargeMpos;
