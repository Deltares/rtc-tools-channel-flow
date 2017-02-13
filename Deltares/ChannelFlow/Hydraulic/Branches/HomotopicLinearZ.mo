within Deltares.ChannelFlow.Hydraulic.Branches;

model HomotopicLinearZ
  import SI = Modelica.SIunits; 
  extends Deltares.ChannelFlow.Hydraulic.Branches.HomotopicLinear(redeclare connector HQPort1 = Deltares.ChannelFlow.Interfaces.HQZCPort);  
  SI.VolumeFlowRate[n_level_nodes+1] Z;
  SI.VolumeFlowRate[n_level_nodes] C;
equation

  for section in 1:n_level_nodes loop
    theta*_cross_section[section] * der(C[section])+(1-theta)*(nominal_width[section] * nominal_depth[section] *der(C[section]) ) 
    -(Z[section+1]  - Z[section])/dx=0;    
  end for;

  
//  if(Q[1] > 0)then
   for section in 2:n_level_nodes loop
      Z[section] = theta*Q[section] * C[section-1]+(1-theta) *(Q[section] * 1.0 +Q_nominal*C[section-1] - Q_nominal*1.0 );
    end for;

//  else
//     for section in 2:n_level_nodes loop
//        Z[section] = theta*Q[section] * C[section]+(1-theta) *(Q[section] * 1.0 +Q_nominal*C[section] - Q_nominal*1.0 );        
//    end for;

//  end if;
    Z[1] = HQUp.Z;
    Z[n_level_nodes+1] = -HQDown.Z;   
    HQDown.C  = C[n_level_nodes] ;
    HQUp.C    = C[1] ;
  
end HomotopicLinearZ;