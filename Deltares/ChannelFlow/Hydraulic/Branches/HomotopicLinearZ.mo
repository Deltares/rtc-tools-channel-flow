within Deltares.ChannelFlow.Hydraulic.Branches;

model HomotopicLinearZ
  import SI = Modelica.SIunits; 
  extends Deltares.ChannelFlow.Hydraulic.Branches.HomotopicLinear(redeclare connector HQPort = Deltares.ChannelFlow.Interfaces.HQCMPort);
  SI.VolumeFlowRate[n_level_nodes + 1] M;
  SI.Density[n_level_nodes] C(each min = 0);
  parameter Real C_nominal = 1;
  parameter SI.Distance dx2 = length / (n_level_nodes);
equation
  //A*dc/dt + dz/dx=0, this equation is linearized
  for section in 1:n_level_nodes loop
    theta * der(_cross_section[section] * C[section]) + (1 - theta) * (nominal_width[section] * nominal_depth[section] * der(C[section]) + C_nominal * der(_cross_section[section])) + (M[section + 1] - M[section]) / dx2 = 0;    
  end for; 
  //calculation of the salt mass flow rate at the internal boundary points
  for section in 2:n_level_nodes loop
    if(Q[section] > 0)then
      M[section] = theta * Q[section] * C[section - 1] + (1 - theta) *(Q_nominal * C_nominal + C_nominal * (Q[section] - Q_nominal) + Q_nominal * (C[section - 1] - C_nominal));
    elseif (Q[section] < 0)then
      M[section] = theta * Q[section] * C[section] + (1 - theta) * (-Q_nominal * C_nominal + C_nominal * (Q[section] + Q_nominal) - Q_nominal * (C[section] - C_nominal));
    else
      M[section] = 0;
    end if;
  end for;
  //setting of the salt mass flow rate and the concentration of salt at the connections.
  M[1] = HQUp.M;
  M[n_level_nodes + 1] = -HQDown.M;
  HQDown.C = C[n_level_nodes];
  HQUp.C =C[1];
end HomotopicLinearZ;
