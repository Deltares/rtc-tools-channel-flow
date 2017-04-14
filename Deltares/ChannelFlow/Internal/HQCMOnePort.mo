within Deltares.ChannelFlow.Internal;

partial model HQCMOnePort "Partial model of one port"
  replaceable package medium = Deltares.media.FreshWater;
  Deltares.ChannelFlow.Interfaces.HQCMPort HQCM(redeclare package medium = medium) annotation(Placement(visible = true, transformation(origin = {0, -80}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {0, -80}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
 // Deltares.ChannelFlow.Interfaces.HQCMPort HQ = HQCM; // For backwards compatibility
end HQCMOnePort;
