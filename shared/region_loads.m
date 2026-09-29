function loads = region_loads(lesions, region_masks)
%REGION_LOADS Proportion of each region destroyed, for every patient.
%   loads = region_loads(lesions, region_masks)
%   lesions: N x V logical; region_masks: R x V logical (same voxels).
%   loads: N x R, 0-1. Empty regions give 0.
sizes = sum(region_masks, 2)';
loads = (double(lesions) * double(region_masks')) ./ max(sizes, 1);
end
