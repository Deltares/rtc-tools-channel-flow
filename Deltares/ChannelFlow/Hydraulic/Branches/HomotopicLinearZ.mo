within Deltares.ChannelFlow.Hydraulic.Branches;

model HomotopicLinearZ
  import SI = Modelica.SIunits; 
  extends Deltares.ChannelFlow.Hydraulic.Branches.HomotopicLinear(redeclare connector HQPort = Deltares.ChannelFlow.Interfaces.HQZCPort);  
  SI.VolumeFlowRate[n_level_nodes+1] Z;
  SI.Density[n_level_nodes] C;
  parameter Real C_nominal=1;
equation

  for section in 1:n_level_nodes loop
    theta*_cross_section[section] * der(C[section])+(1-theta)*(nominal_width[section] * nominal_depth[section] *der(C[section]) ) 
    +(Z[section+1]  - Z[section])/dx=0;    
  end for;  

   for section in 2:n_level_nodes loop
      if(Q[section] > 0)then
            Z[section] = (theta*Q[section] * C[section-1]+(1-theta) *(Q_nominal*C_nominal + C_nominal*(Q[section]-Q_nominal) +              Q_nominal*(C[section-1]-C_nominal)));
      else 
      Z[section] = 
      (theta*Q[section] * C[section]+(1-theta) *(-Q_nominal*C_nominal + C_nominal*(Q[section]+Q_nominal) - Q_nominal*(C[section]-C_nominal)));
       end if;
    end for;

    Z[1] = HQUp.Z;
    Z[n_level_nodes+1] = -HQDown.Z;   
    HQDown.C  = C[n_level_nodes] ;
    HQUp.C    = C[1] ;    

end HomotopicLinearZ;
