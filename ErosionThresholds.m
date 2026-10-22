function [Ecr, Ecs] = ErosionThresholds(p)

% [Ecr, Ecs] = ErosionThresholds(p)
%
% Threshold erosion rates (m/yr) of the smoothed erosion law
% section 14) for bedrock (Ecr) and sediment (Ecs): p.Ecr and p.Ecs with
% p.channelThreshold = 'smooth', zero with 'area' (no threshold in the law;
% erosion only where A > p.Ac) and Ecs = 0 without sediment.

Ecr = 0; Ecs = 0;
if strcmp(p.channelThreshold, 'smooth')
    Ecr = p.Ecr;
    if p.doSediment, Ecs = p.Ecs; end
end

end
