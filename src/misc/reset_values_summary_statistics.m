%% assign all the values to each vector


spine_dendrite_type_label5 = layer5_dendrite_level1_summary.spine_dendrite_type_label; % [1: filopodia, 2: mushroom, 3: long thin, 4: thin, 5: stubby, dendrite_label] Ns x 2
dendrite_length5  = layer5_dendrite_level1_summary.dendrite_length; % Nd x 1
dendrite_radius5 = layer5_dendrite_level1_summary.dendrite_radius; % Nd x 1 
dendrite_each_spine_type_number5 = layer5_dendrite_level1_summary.dendrite_each_spine_type_number; % specifically add the 6th type as the branched spine, simply to check between different layers Nd x 6
neuron_total_type_number5 = layer5_dendrite_level1_summary.neuron_total_type_number; % Nn x 6

spine_dendrite_type_label4 = layer4_dendrite_level1_summary.spine_dendrite_type_label; % [1: filopodia, 2: mushroom, 3: long thin, 4: thin, 5: stubby, dendrite_label] Ns x 2
dendrite_length4  = layer4_dendrite_level1_summary.dendrite_length; % Nd x 1
dendrite_radius4 = layer4_dendrite_level1_summary.dendrite_radius; % Nd x 1 
dendrite_each_spine_type_number4 = layer4_dendrite_level1_summary.dendrite_each_spine_type_number; % specifically add the 6th type as the branched spine, simply to check between different layers Nd x 6
neuron_total_type_number4 = layer4_dendrite_level1_summary.neuron_total_type_number; % Nn x 6

spine_dendrite_type_label2 = layer2_dendrite_level1_summary.spine_dendrite_type_label; % [1: filopodia, 2: mushroom, 3: long thin, 4: thin, 5: stubby, dendrite_label] Ns x 2
dendrite_length2  = layer2_dendrite_level1_summary.dendrite_length; % Nd x 1
dendrite_radius2 = layer2_dendrite_level1_summary.dendrite_radius; % Nd x 1 
dendrite_each_spine_type_number2 = layer2_dendrite_level1_summary.dendrite_each_spine_type_number; % specifically add the 6th type as the branched spine, simply to check between different layers Nd x 6
neuron_total_type_number2 = layer2_dendrite_level1_summary.neuron_total_type_number; % Nn x 6




singleSynHeadVolume5 = layer5_summary.singleSynHeadVolume;
singleSynMeanHeadRadius5 = layer5_summary.singleSynMeanHeadRadius;
singleSynNeckLength5 = layer5_summary.singleSynNeckLength;
singleSynNeckSection5 = layer5_summary.singleSynNeckSection;
singleSynNeckMeanRadius5 = layer5_summary.singleSynNeckMeanRadius;
doubleSynHeadVolume5 = layer5_summary.doubleSynHeadVolume;
doubleSynMeanHeadRadius5 = layer5_summary.doubleSynMeanHeadRadius;
doubleSynNeckLength5 = layer5_summary.doubleSynNeckLength;
doubleSynNeckSection5 = layer5_summary.doubleSynNeckSection;
doubleSynNeckMeanRadius5 = layer5_summary.doubleSynNeckMeanRadius;
singleSynapticCleftSize5 = layer5_summary.singleSynapticCleftSize;
% singleDendriteDensity5 = layer5_summary.singleDendriteDensity;
sinsperimeterRatio5 = layer5_summary.sinsperimeterRatio;
sinsperimeterWeightedWrappingArea5 = layer5_summary.sinsperimeterWeightedWrappingArea;
sinspostSynapseTouchingArea5 = layer5_summary.sinspostSynapseTouchingArea;
sinspostSynapseTouchingRatio5 = layer5_summary.sinspostSynapseTouchingRatio;
sinspreSynapseTouchingArea5 = layer5_summary.sinspreSynapseTouchingArea;
sinspreSynapseTouchingRatio5 = layer5_summary.sinspreSynapseTouchingRatio;
singleSynHeadNeckTouchingArea5 = layer5_summary.singleSynHeadNeckTouchingArea;
singleSynHeadNeckTouchingRatio5 = layer5_summary.singleSynHeadNeckTouchingRatio;
% doubleDendriteDensity5 = layer5_summary.doubleDendriteDensity;
doubleSynapticCleftSize5 = layer5_summary.doubleSynapticCleftSize;
dousperimeterRatio5 = layer5_summary.dousperimeterRatio;
dousperimeterWeightedWrappingArea5 = layer5_summary.dousperimeterWeightedWrappingArea;
douspostSynapseTouchingArea5  = layer5_summary.douspostSynapseTouchingArea;
douspostSynapseTouchingRatio5 = layer5_summary.douspostSynapseTouchingRatio;
douspreSynapseTouchingArea5  = layer5_summary.douspreSynapseTouchingArea;
douspreSynapseTouchingRatio5 = layer5_summary.douspreSynapseTouchingRatio;
doubleSynHeadNeckTouchingArea5 = layer5_summary.doubleSynHeadNeckTouchingArea;
doubleSynHeadNeckTouchingRatio5 = layer5_summary.doubleSynHeadNeckTouchingRatio;
% doubleSynPosition5 = layer5_summary.doubleSynPosition;
% spineCoordinate5  =layer5_summary.spineCoordinate;





singleSynHeadVolume2 = layer2_summary.singleSynHeadVolume;
singleSynMeanHeadRadius2 = layer2_summary.singleSynMeanHeadRadius;
singleSynNeckLength2 = layer2_summary.singleSynNeckLength;
singleSynNeckSection2 = layer2_summary.singleSynNeckSection;
singleSynNeckMeanRadius2 = layer2_summary.singleSynNeckMeanRadius;
doubleSynHeadVolume2 = layer2_summary.doubleSynHeadVolume;
doubleSynMeanHeadRadius2 = layer2_summary.doubleSynMeanHeadRadius;
doubleSynNeckLength2 = layer2_summary.doubleSynNeckLength;
doubleSynNeckSection2 = layer2_summary.doubleSynNeckSection;
doubleSynNeckMeanRadius2 = layer2_summary.doubleSynNeckMeanRadius;
singleSynapticCleftSize2 = layer2_summary.singleSynapticCleftSize;
sinsperimeterRatio2 = layer2_summary.sinsperimeterRatio;
sinsperimeterWeightedWrappingArea2 = layer2_summary.sinsperimeterWeightedWrappingArea;
sinspostSynapseTouchingArea2 = layer2_summary.sinspostSynapseTouchingArea;
sinspostSynapseTouchingRatio2 = layer2_summary.sinspostSynapseTouchingRatio;
sinspreSynapseTouchingArea2 = layer2_summary.sinspreSynapseTouchingArea;
sinspreSynapseTouchingRatio2 = layer2_summary.sinspreSynapseTouchingRatio;
singleSynHeadNeckTouchingArea2 = layer2_summary.singleSynHeadNeckTouchingArea;
singleSynHeadNeckTouchingRatio2 = layer2_summary.singleSynHeadNeckTouchingRatio;
doubleSynapticCleftSize2 = layer2_summary.doubleSynapticCleftSize;
dousperimeterRatio2 = layer2_summary.dousperimeterRatio;
dousperimeterWeightedWrappingArea2 = layer2_summary.dousperimeterWeightedWrappingArea;
douspostSynapseTouchingArea2  = layer2_summary.douspostSynapseTouchingArea;
douspostSynapseTouchingRatio2 = layer2_summary.douspostSynapseTouchingRatio;
douspreSynapseTouchingArea2 = layer2_summary.douspreSynapseTouchingArea;
douspreSynapseTouchingRatio2 = layer2_summary.douspreSynapseTouchingRatio;
doubleSynHeadNeckTouchingArea2 = layer2_summary.doubleSynHeadNeckTouchingArea;
doubleSynHeadNeckTouchingRatio2 = layer2_summary.doubleSynHeadNeckTouchingRatio;
% doubleSynPosition2 = layer2_summary.doubleSynPosition;
% spineCoordinate2  =layer2_summary.spineCoordinate;
% singleDendriteDensity2 = layer2_summary.singleDendriteDensity;
% doubleDendriteDensity2 = layer2_summary.doubleDendriteDensity;

singleSynHeadVolume4 = layer4_summary.singleSynHeadVolume;
singleSynMeanHeadRadius4 = layer4_summary.singleSynMeanHeadRadius;
singleSynNeckLength4 = layer4_summary.singleSynNeckLength;
singleSynNeckSection4 = layer4_summary.singleSynNeckSection;
singleSynNeckMeanRadius4 = layer4_summary.singleSynNeckMeanRadius;
doubleSynHeadVolume4 = layer4_summary.doubleSynHeadVolume;
doubleSynMeanHeadRadius4 = layer4_summary.doubleSynMeanHeadRadius;
doubleSynNeckLength4 = layer4_summary.doubleSynNeckLength;
doubleSynNeckSection4 = layer4_summary.doubleSynNeckSection;
doubleSynNeckMeanRadius4 = layer4_summary.doubleSynNeckMeanRadius;
singleSynapticCleftSize4 = layer4_summary.singleSynapticCleftSize;
sinsperimeterRatio4 = layer4_summary.sinsperimeterRatio;
sinsperimeterWeightedWrappingArea4 = layer4_summary.sinsperimeterWeightedWrappingArea;
sinspostSynapseTouchingArea4 = layer4_summary.sinspostSynapseTouchingArea;
sinspostSynapseTouchingRatio4 = layer4_summary.sinspostSynapseTouchingRatio;
sinspreSynapseTouchingArea4 = layer4_summary.sinspreSynapseTouchingArea;
sinspreSynapseTouchingRatio4 = layer4_summary.sinspreSynapseTouchingRatio;
singleSynHeadNeckTouchingArea4 = layer4_summary.singleSynHeadNeckTouchingArea;
singleSynHeadNeckTouchingRatio4 = layer4_summary.singleSynHeadNeckTouchingRatio;
doubleSynapticCleftSize4 = layer4_summary.doubleSynapticCleftSize;
dousperimeterRatio4 = layer4_summary.dousperimeterRatio;
dousperimeterWeightedWrappingArea4 = layer4_summary.dousperimeterWeightedWrappingArea;
douspostSynapseTouchingArea4  = layer4_summary.douspostSynapseTouchingArea;
douspostSynapseTouchingRatio4 = layer4_summary.douspostSynapseTouchingRatio;
douspreSynapseTouchingArea4  = layer4_summary.douspreSynapseTouchingArea;
douspreSynapseTouchingRatio4 = layer4_summary.douspreSynapseTouchingRatio;
doubleSynHeadNeckTouchingArea4 = layer4_summary.doubleSynHeadNeckTouchingArea;
doubleSynHeadNeckTouchingRatio4 = layer4_summary.doubleSynHeadNeckTouchingRatio;