% Quantification of the dendrite branches and spines for the different layers based on the dataset extracted from MICrONS
function main_plot_quantification(data_path, save_path)

% 5. Summarize the results (plots...)

% % example plot the quantifications 
multilayer_summary_path = fullfile(data_path, 'layer2345_quantification_result_20240813'); % path to the summary of the quantification for the different layers
% % read in the summary from the quantification for each layer
% % both the quantification of the dendrite branches as well the quantification for the dendrite spines.
plot_save_folder = save_path; % path to the folder where the plots should be saved
[layer23_dendrite_summary, layer23_spine_summary, layer4_dendrite_summary, layer4_spine_summary, layer5_dendrite_summary, layer5_spine_summary] = structQuant.read_in_summary(multilayer_summary_path);
structQuant.plot_summary_plot_example(save_path,layer23_dendrite_summary, layer23_spine_summary, layer4_dendrite_summary, layer4_spine_summary, layer5_dendrite_summary, layer5_spine_summary);

