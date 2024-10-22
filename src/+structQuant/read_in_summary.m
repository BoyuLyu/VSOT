function [layer23_dendrite_summary, layer23_spine_summary, layer4_dendrite_summary, layer4_spine_summary, layer5_dendrite_summary, layer5_spine_summary] = read_in_summary(multilayer_summary_path)

    layer23_dendrite_summary = load(fullfile(multilayer_summary_path, 'layer2_dendrite_level1_summary.mat'));
    layer23_dendrite_summary = layer23_dendrite_summary.layer2_dendrite_level1_summary;
    layer23_spine_summary = load(fullfile(multilayer_summary_path, 'layer2_summary.mat'));
    layer23_spine_summary = layer23_spine_summary.layer2_summary;
    layer4_dendrite_summary = load(fullfile(multilayer_summary_path, 'layer4_dendrite_level1_summary.mat'));
    layer4_dendrite_summary = layer4_dendrite_summary.layer4_dendrite_level1_summary;
    layer4_spine_summary = load(fullfile(multilayer_summary_path, 'layer4_summary.mat'));
    layer4_spine_summary = layer4_spine_summary.layer4_summary;
    layer5_dendrite_summary = load(fullfile(multilayer_summary_path, 'layer5_dendrite_level1_summary.mat'));
    layer5_dendrite_summary = layer5_dendrite_summary.layer5_dendrite_level1_summary;
    layer5_spine_summary = load(fullfile(multilayer_summary_path, 'layer5_summary.mat'));
    layer5_spine_summary = layer5_spine_summary.layer5_summary;








end