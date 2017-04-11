within Deltares.ChannelFlow.Hydraulic.Branches;

model HomotopicLinear
  import SI = Modelica.SIunits;
  extends Internal.PartialHomotopic(nominal_depth = fill(uniform_nominal_depth, n_level_nodes + 1), nominal_width = linspace(width_up, width_down, n_level_nodes + 1), H(min = H_b));
  // Nominal depth
  parameter SI.Distance uniform_nominal_depth;
  // Upstream Width (same 'Up' as HQUp)
  parameter SI.Distance width_up; 
  // Downstream Width (same 'Down' as HQDown)
  parameter SI.Distance width_down;
  // Array of Widths
  parameter SI.Distance width[n_level_nodes] = linspace(width_up, width_down, n_level_nodes);
  // Upstream Bottom Level (same 'Up' as HQUp)
  parameter SI.Position H_b_up; 
  // Downstream Bottom Level (same 'Down' as HQDown)
  parameter SI.Position H_b_down;
  // Array of Bottom Levels
  parameter SI.Position H_b[n_level_nodes] = linspace(H_b_up, H_b_down, n_level_nodes);
  
  SI.VolumeFlowRate[n_level_nodes + 1,HQCMUp.NOS] M;
  SI.Density[n_level_nodes,HQCMUp.NOS] C(each min = 0);
  parameter Real C_nominal[HQCMUp.NOS] = fill(1,HQCMUp.NOS);
  parameter SI.Distance dx2 = length / (n_level_nodes);
equation
  // Compute cross sections
  _cross_section = width .* (H .- H_b);
  
  //substance transport part
    //A*dc/dt + dz/dx=0, this equation is linearized
  for section in 1:n_level_nodes loop
    theta * der(_cross_section[section] * C[section,:]) = -(1 - theta) * (nominal_width[section] * nominal_depth[section] * der(C[section,:]) + C_nominal * der(_cross_section[section]))- (M[section + 1,:] - M[section,:]) / dx2 ;    
  end for; 
  //calculation of the salt mass flow rate at the internal boundary points
  for section in 2:n_level_nodes loop
    if(Q[section] > 0)then
      M[section,:] = theta .* Q[section] .* C[section - 1,:] + (1 - theta) *(Q_nominal * C_nominal + C_nominal * (Q[section] - Q_nominal) + Q_nominal * (C[section - 1,:] - C_nominal));
    else
      M[section,:] = theta * Q[section] .* C[section,:] + (1 - theta) * (Q_nominal * C_nominal + C_nominal * (Q[section] - Q_nominal) + Q_nominal * (C[section,:] - C_nominal));
    end if;
  end for;
  //setting of the salt mass flow rate and the concentration of salt at the connections.
  M[1,:] = HQCMUp.M;
  M[n_level_nodes + 1,:] = -HQCMDown.M;
  HQCMDown.C = C[n_level_nodes,:];
  HQCMUp.C =C[1,:];
end HomotopicLinear;
