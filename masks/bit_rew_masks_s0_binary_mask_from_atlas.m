%% bit_rew_s0_binary_mask_from_atlas.m
%
% This script creates a binary mask by combining regions from one or more
% atlases, and automatically saves the fmri_mask_image, and .nii versions
% of it in maskdir of your dataset
% 
% USAGE
%
% Script should be run from the root directory of the superdataset, e.g.
% /data/proj_discoverie
%
% DEPENDENCIES
%
% CANlab's CanlabCore and Neuroimaging_Pattern_Masks Github repos on your Matlab path
% if needed, clone from https://github.com/canlab
%
%__________________________________________________________________________
%
% authors: Aleksandra Budzinska, Lukas Van Oudenhove
% date:   KU Leuven, July, 2022
%
%__________________________________________________________________________
% @(#)% LaBGAScore_atlas_binary_mask_from_atlas.m         v1.2       
% last modified: 2022/07/26


% Define maskname and directory where mask will be written
%--------------------------------------------------------------------------

maskname = 'bit_rew_mask_all_regions';
maskdir = '/data/proj_bitter-reward/secondlevel/model_1m_food_images/masks';

% NOTE: in default LaBGAS file organization this is a model-specific dir in
% rootdir/secondlevel/model_xxx, but you can obviously specify any dir here


% Load atlas of your choice
%--------------------------------------------------------------------------
canlab = load_atlas('canlab2023_coarse_2mm'); % you get a CANlab atlas object

% NOTE: you can not only load CANlab atlases using keywords as in this example,
% but any atlas you like by creating a new atlas object from a .nii atlas image, 
% e.g. aicha = atlas('/opt/KUL_apps/spm12/atlas/AICHA.nii');
% see help atlas for adding labels etc


% Select regions you want to include from the selected atlas 
%--------------------------------------------------------------------------
canlab_subset = select_atlas_subset(canlab, {'Ctx_p24pr_L','Ctx_p24pr_R','Ctx_33pr_L','Ctx_33pr_R','Ctx_a24pr_L','Ctx_a24pr_R','Ctx_p32pr_L','Ctx_p32pr_R','Ctx_a24_L','Ctx_a24_R','Ctx_d32_L','Ctx_d32_R','Ctx_8BM_L','Ctx_8BM_R','Ctx_p32_L','Ctx_p32_R','Ctx_10r_L','Ctx_10r_R','Ctx_47m_L','Ctx_47m_R','Ctx_10d_L','Ctx_10d_R','Ctx_a47r_L','Ctx_a47r_R','Ctx_a10p_L','Ctx_a10p_R','Ctx_10pp_L','Ctx_10pp_R','Ctx_11l_L','Ctx_11l_R','Ctx_13l_L','Ctx_13l_R','Ctx_OFC_L','Ctx_OFC_R','Ctx_47s_L','Ctx_47s_R','Ctx_52_L','Ctx_52_R','Ctx_PFcm_L','Ctx_PFcm_R','Ctx_PoI2_L','Ctx_PoI2_R','Ctx_FOP4_L','Ctx_FOP4_R','Ctx_MI_L','Ctx_MI_R','Ctx_Pir_L','Ctx_Pir_R','Ctx_AVI_L','Ctx_AVI_R','Ctx_AAIC_L','Ctx_AAIC_R','Ctx_FOP1_L','Ctx_FOP1_R','Ctx_FOP3_L','Ctx_FOP3_R','Ctx_FOP2_L','Ctx_FOP2_R','Ctx_25_L','Ctx_25_R','Ctx_s32_L','Ctx_s32_R','Ctx_pOFC_L','Ctx_pOFC_R','Ctx_PoI1_L','Ctx_PoI1_R','Ctx_Ig_L','Ctx_Ig_R','Ctx_FOP5_L','Ctx_FOP5_R','Ctx_p10p_L','Ctx_p10p_R','Ctx_PI_L','Ctx_PI_R','Ctx_a32pr_L','Ctx_a32pr_R','Ctx_p24_L','Ctx_p24_R',...
    'NAc_R','NAc_L','CAU_R','CAU_L','PUT_R','PUT_L','GP_R','GP_L',...
    'Hythal_L','Hythal_R','Haben_L','Haben_R','Subic_L','Hippocampus_L','Amygdala_L','Subic_R','Hippocampus_R','Amygdala_R','BST_SLEA_L','BST_SLEA_R',...
    'Bstem_PAG','Bstem_SN_L','Bstem_SN_R','Bstem_VTA_L','Bstem_VTA_L','Bstem_PBP_L','Bstem_PBP_L','Bstem_VSM_L','Bstem_VSM_R','Bstem_RN_L','Bstem_RN_R','Bstem_LCSubC_L','Bstem_LCSubC_R','Bstem_PBN_L','Bstem_PBN_R'});

% NOTE: you can also use the numbers rather than the labels of regions to
% select as a vector
% see help atlas.select_atlas_subset

canlab_subset = canlab_subset.threshold(0.20); % canlab2023 is a probabilistic atlas


% Load and select the regions from another atlas (e.g. subcortical) like above
%--------------------------------------------------------------------------
% subcortical = load_atlas('cit168'); % you get a CANlab atlas object
% subcortical_subset = select_atlas_subset(subcortical, {'Put','Cau','NAC','BST_SLEA','GPe','GPi','SNc','RN','SNr','PBP','VTA','VeP','Hythal','STN'}); % still a CANlab atlas object


% Change data type of probability maps if needed
%--------------------------------------------------------------------------

% NOTE: Before merging the atlases, make sure that they have the same type of variables (look at the probability_maps. E.g. single/sparse double)
% This command converts probability_maps property variable type from single to double, and then to sparse double - may not be needed for all atlases

% subcortical_subset.probability_maps = double(subcortical_subset.probability_maps);
% subcortical_subset.probability_maps = sparse(subcortical_subset.probability_maps);


% Merge the atlases
%--------------------------------------------------------------------------

combined_atlas = canlab_subset; 

% NOTE: merge_atlases function automatically resamples the space of the second to the first, and adds consecutive labels


% Convert the atlas object to a CANlab fmri_mask_image object
%--------------------------------------------------------------------------

mask = fmri_mask_image(canlab_subset);

% NOTE: this automatically binarizes the mask!


% Save your atlas subset & binarized mask in maskdir
%--------------------------------------------------------------------------

write(mask, 'fname', fullfile(maskdir,[maskname,'.nii'])); % writes mask to .nii file
write(combined_atlas, 'fname', fullfile(maskdir,[maskname,'_atlas.nii'])); % writes mask to .nii file

savefilenamedata = fullfile(maskdir,[maskname,'.mat']); % saves fmri_mask_image object as .mat file
save(savefilenamedata, 'mask', 'combined_atlas','-v7.3');
fprintf('\nSaved mask\n');

