within Deltares.ChannelFlow.Hydraulic.Reservoir;


model HomotopicVolume
  // This model introduces a nonlinear, non-convex relationship between
  // reservoir volume and water level using a fourth-order polynomial.
  //
  // The optimisation starts from a linear approximation,
  // water_level = volume / A, where A is the reservoir surface area.
  //
  // During homotopy optimisation, the parameter `theta` gradually
  // transitions the model from the linear approximation to the full
  // nonlinear polynomial relation.
  //
  // This model must be used together with the HomotopyMixin.

  extends Internal.PartialHomotopicVolume;

end HomotopicVolume;

