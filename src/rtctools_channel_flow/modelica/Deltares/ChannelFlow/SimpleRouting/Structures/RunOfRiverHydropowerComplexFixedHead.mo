within Deltares.ChannelFlow.SimpleRouting.Structures;

block RunOfRiverHydropowerComplexFixedHead "Node for a simple complex of a run-of-river hydropower plant and a weir. Head difference for power production is constant."
  import SI = Modelica.SIunits;
  extends Deltares.ChannelFlow.Internal.QSISO;
    // Head difference
  parameter SI.Position dH;
  // Turbine efficiency
  parameter real nu;
  // Water density
  parameter SI.Density ro;
  // Turbine flow
  output SI.VolumeFlowRate Q_turbine(min=0);
  // Spill flow
  output SI.VolumeFlowRate Q_spill(min=0);
  // Power production
  output SI.Power P;
  equation
    QOut.Q = Q_turbine + Q_spill;
    QOut.Q = QIn.Q;
    P = nu * ro * Deltares.Constants.g_n * dH * Q_turbine;
end RunOfRiverHydropowerComplexFixedHead;
