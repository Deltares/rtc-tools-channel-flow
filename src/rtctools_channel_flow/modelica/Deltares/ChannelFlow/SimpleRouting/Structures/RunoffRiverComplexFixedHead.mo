within Deltares.ChannelFlow.SimpleRouting.Structures;

block RunoffRiverComplexFixedHead "Node for a simple complex of a runoff-river power plant and a weir. Head difference for power production is constant."
  extends Deltares.ChannelFlow.Internal.QSISO;
  import SI = Modelica.SIunits;
  input SI.VolumeFlowRate Q_turbine;
  output SI.VolumeFlowRate Q_spill(min=0);
  Real dH;
  Real efficiency;
  Real Power;
  Real density;
  equation
    QOut.Q = Q_turbine + Q_spill;
    QOut.Q = QIn.Q;
    Power = Q_turbine * efficiency * density * dH;
end RunoffRiverComplexFixedHead;
