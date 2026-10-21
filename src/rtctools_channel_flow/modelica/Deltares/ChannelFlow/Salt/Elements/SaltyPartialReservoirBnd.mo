within Deltares.ChannelFlow.Salt.Elements;

partial model SaltyPartialReservoirBnd
  extends Deltares.ChannelFlow.Internal.HQTwoPort;
  extends Deltares.ChannelFlow.Internal.QForcing;
  extends Deltares.ChannelFlow.Internal.QLateral;
  extends Deltares.ChannelFlow.Internal.Reservoir;
  // States
  Modelica.Units.SI.Position H;
  Modelica.Units.SI.Concentration C(nominal = 0.001);

  parameter Real theta = 1.0;
  parameter Modelica.Units.SI.VolumeFlowRate Q_nominal = 1.0;
  parameter Modelica.Units.SI.Concentration C_nominal = 1.0;
  parameter Modelica.Units.SI.Volume V_nominal;

equation
  // Water level
  H = HQUp.H;
  H = HQDown.H;

  C = HQUp.C[1];
  C = HQDown.C[1];

  // Mass balance
  der(V) = 0; 
  der(V * C) = 0;


end SaltyPartialReservoirBnd;
