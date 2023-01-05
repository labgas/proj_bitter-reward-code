% LCN12_JULIE_first_level_analysis.m
%
% This script will perform the first level analysis for the bitter-reward 
% study of Julie. We assume that the data are preprocessed and organized
% in the following way:
%   Each subject has his own folder. Within this folder, there should be
%       a subfolder processed which contains the preprocessed data
%                             (smoothed and warped 4D fMRI images) according 
%                             to the naming convention s8waXXX_func_YYY_runZ.nii 
%                             in which XXX is the subjectcode, YYY is the 
%                             condition (con or bit) and Z is the run 
%                             number (1-4). It also contains the
%                             realignment files (one file per run) with the
%                             naming rp_XXX_func_YYY_runZ.txt
%       a subfolder logfiles which contains for each run the logfiles to
%                            extract the timings. The naming convention is 
%                            as follows: XXX_YYY_runZ_photo_results.txt for 
%                            run 1 and 2 and XXX_YYY_runZ-FID_1_*.log for 
%                            run 3 and 4.  
% A subfolder results will be created in which you find the first level 
% analysis for each of the two conditions. The batch files will be written 
% in the subfolder batch (found in the main directory of each subject).
%
% The design of the experiment is as follows:
%   There are two sessions (between brackets the coding): 
%         bitter (bit) or placebo (con).
%   Each session had 2 x 2 runs in which the first two runs consist of
%         passive viewing of food images and the last two runs are a food
%         incentive delay task. 
%
%         The passive viewing of food images consist of a block design in 
%         which there are 5 conditions: 
%         1) a rest block of 24s, 
%         2-4) a block with 3 possible sorts of images:  high calory
%         images (HC), low calory images (LC) and non-food images (NF). A
%         block of images consist of the presentation of 8 different images
%         for 3s so these blocks also last for 24s. 
%         5) After each image block a rating block was present (i.e. the 
%         fifth condition). 
%         The timings can be derived from the logfile XXX_YYY_runZ_photo_results.txt.
%         Important to note is that run 1 contains 216 volumes and run 2
%         contains 228 volumes. However, in the last part of run2 a GI
%         rating was performed and we don't need this in our current
%         analysis. Therefore, we neglect the last 12 volumes.
%
%         The food incentive delay task (FID) consist of a cue (3
%         possible cues C0, C2 and C10) for 750ms followed by a delay of
%         3s. Then a target was shown for 1s (arrow to the left or to the right)
%         and the subject has to press one of two buttons corresponding to 
%         each side. The button press had to be within the interval that 
%         the target was on the screen (i.e. within 1s after target onset).
%         If pressed correctly, the subject gained 0, 2 or 10 "points" 
%         depending on the initial cue C0, C2 or C10. Feedback was shown 
%         for 1.5s on the screen in the form of a picture. If pressed 
%         incorrect, the subject would gain nothing and if pressed too late
%         or when omitted, the subject would lose 10 "points" irrespective 
%         of the initial cue. This latter condition as well as the missed 
%         events were rare and will not be modelled. However, in 15 of the 
%         55 trials, the subjects would not gain anything despite the fact 
%         that they answered correctly. These events will be modelled as 
%         as well. We will model the following conditions as events:
%         1) cue C0
%         2) cue C2
%         3) cue C10
%         4) feedback C10 win
%         5) feedback C10 nowin
%         6) feedback C2 win
%         7) feedback C2 nowin
%         8) feedback C0
%
% Each subject always got for the FID task the same order of conditions
% based on the presentation scripts. There were 55 trials in total divided
% as follows:
%   12 trials feedback C10 win
%   8  trials feedback C10 nowin
%   13 trials feedback C2 win
%   7 trials feedback C2 nowin
%   15 trials feedback C0
%
% Note: trials in which the subject responded wrong, too late or omitted to
% respond, were not modelled (too few instances).
%
% author:   Patrick Dupont
% date:     June 2018
% history:  July 2018: presentation logfile was not always consistent at 
%                      the last columns. Only the first 5 columns are 
%                      actually used in the script so the rest will no 
%                      longer read + bugfix when run 3 was missing but not
%                      run 4.rre
%__________________________________________________________________________
% LCN12_JULIE_first_level_analysis.m      v0.11   last modified: 2018/07/01

clear
close all

SUBJECTLIST = {
    'J:\GBW-0264_TARGID-Brain-Gut-Axis\JULIE\Bitter_reward_Julie\Brain_responses\subj01'
    'J:\GBW-0264_TARGID-Brain-Gut-Axis\JULIE\Bitter_reward_Julie\Brain_responses\subj02'
    'J:\GBW-0264_TARGID-Brain-Gut-Axis\JULIE\Bitter_reward_Julie\Brain_responses\subj04'
    'J:\GBW-0264_TARGID-Brain-Gut-Axis\JULIE\Bitter_reward_Julie\Brain_responses\subj05'
    'J:\GBW-0264_TARGID-Brain-Gut-Axis\JULIE\Bitter_reward_Julie\Brain_responses\subj06'
    'J:\GBW-0264_TARGID-Brain-Gut-Axis\JULIE\Bitter_reward_Julie\Brain_responses\subj07'
    'J:\GBW-0264_TARGID-Brain-Gut-Axis\JULIE\Bitter_reward_Julie\Brain_responses\subj10'
    'J:\GBW-0264_TARGID-Brain-Gut-Axis\JULIE\Bitter_reward_Julie\Brain_responses\subj08'
    'J:\GBW-0264_TARGID-Brain-Gut-Axis\JULIE\Bitter_reward_Julie\Brain_responses\subj09'
    'J:\GBW-0264_TARGID-Brain-Gut-Axis\JULIE\Bitter_reward_Julie\Brain_responses\subj11'
    'J:\GBW-0264_TARGID-Brain-Gut-Axis\JULIE\Bitter_reward_Julie\Brain_responses\subj12'
    'J:\GBW-0264_TARGID-Brain-Gut-Axis\JULIE\Bitter_reward_Julie\Brain_responses\subj13'
    'J:\GBW-0264_TARGID-Brain-Gut-Axis\JULIE\Bitter_reward_Julie\Brain_responses\subj14'
    'J:\GBW-0264_TARGID-Brain-Gut-Axis\JULIE\Bitter_reward_Julie\Brain_responses\subj15'
    'J:\GBW-0264_TARGID-Brain-Gut-Axis\JULIE\Bitter_reward_Julie\Brain_responses\subj16'
    'J:\GBW-0264_TARGID-Brain-Gut-Axis\JULIE\Bitter_reward_Julie\Brain_responses\subj17'
    'J:\GBW-0264_TARGID-Brain-Gut-Axis\JULIE\Bitter_reward_Julie\Brain_responses\subj18'
    'J:\GBW-0264_TARGID-Brain-Gut-Axis\JULIE\Bitter_reward_Julie\Brain_responses\subj19'
    'J:\GBW-0264_TARGID-Brain-Gut-Axis\JULIE\Bitter_reward_Julie\Brain_responses\subj20'
    'J:\GBW-0264_TARGID-Brain-Gut-Axis\JULIE\Bitter_reward_Julie\Brain_responses\subj22'
    'J:\GBW-0264_TARGID-Brain-Gut-Axis\JULIE\Bitter_reward_Julie\Brain_responses\subj23'
    };

template_prefix = 's8wa';
template_log1 = '_photo_results.txt';
template_log2 = '-FID_1_final*.log';

nr_runs = 4;
nr_volumes1 = 216; % for run1 and run2
nr_volumes2 = 228; % for run3 and run4

TR = 2.5; % in seconds

% conditions run1 and run2
conditions_name12 = {
    'rest'
    'HC'
    'LC'
    'NF'
    'rating'
    };
% conditions run3 and run4
conditions_name34 = {
    'cue C0'
    'cue C2'
    'cue C10'
    'feedback_C10_win'
    'feedback_C10_nowin'
    'feedback_C2_win'
    'feedback_C2_nowin'
    'feedback_C0'
};

% 1 = feedback_C10_win, 2 = feedback_C10_nowin, 3 = feedback_C2_win
% 4 = feedback_C2_nowin, 5 = feedback_C0
order_trials = [5 4 3 3 1 4 1 5 2 5 ...
                3 3 1 3 2 5 3 2 5 1 ... 
                3 3 1 4 4 5 3 1 3 5 ...
                4 1 5 4 3 5 5 1 1 2 ...
                1 5 2 5 3 2 4 5 5 1 ...
                2 2 5 1 3];
            

%+++++++ END OF SETTINGS ++++++++++++++++++++++++++++++++++++++++++++++++++

nr_subjects = size(SUBJECTLIST,1);
nr_conditions12 = size(conditions_name12,1);
nr_conditions34 = size(conditions_name34,1);

% define settings of SPM
spm('Defaults','fMRI')

for subj = 1:nr_subjects
    clear subjectdir
    
    subjectdir = char(SUBJECTLIST{subj});
    [~,subjectcode,~] = fileparts(subjectdir);
    batchdir = fullfile(subjectdir,'batch');
    datadir  = fullfile(subjectdir,'processed');
    logdir   = fullfile(subjectdir,'logfiles');

    cd(subjectdir)
    mkdir('results')
    cd('results')
    mkdir('bitter')
    mkdir('placebo')
    
    outputdir_bit = fullfile(fullfile(subjectdir,'results'),'bitter');
    outputdir_con = fullfile(fullfile(subjectdir,'results'),'placebo');

    % FOR THE CON STUDY
    %++++++++++++++++++
    
    go1 = 0; % for run1
    go2 = 0; % for run2
    go3 = 0; % for run3
    go4 = 0; % for run4
    
    % FIRST BATCH MODEL SPECIFICATION
    %--------------------------------
    % create empty structures
    matlabbatch = struct([]);
    
    matlabbatch{1}.spm.stats.fmri_spec.dir = {outputdir_con};
    matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs';
    matlabbatch{1}.spm.stats.fmri_spec.timing.RT = TR;
    matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 16;
    matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 1;
    batch_index = 0;
    for sess = 1:2
       % read the scans
        filename = fullfile(datadir,[template_prefix subjectcode '_func_con_run' num2str(sess) '.nii']);
        if exist(filename,'file') == 2
           batch_index = batch_index+1;
           if sess == 1
              go1 = 1;
           elseif sess == 2
              go2 = 1;
           end
           filelist = {};
           for i = 1:nr_volumes1
               filelist{i,1} = fullfile(datadir,[template_prefix subjectcode '_func_con_run' num2str(sess) '.nii,' num2str(i)]);
           end
           matlabbatch{1}.spm.stats.fmri_spec.sess(batch_index).scans = filelist;
             
           % read the logfile
           cd(logdir)
           tmp = dir([subjectcode '_con_run' num2str(sess) '*' template_log1]);        
           fid   = fopen(char(tmp(1).name),'r'); 
           datafile = textscan(fid,'%s %s %s %s %s %s %s %s %s %s %*[^\n]','delimiter','\t');
           fclose(fid);
           % first row of datafile is a header and should be neglected
           % column 8 of datafile is the start time in ms since start of fMRI scan
           % column 9 of datafile is the end time in ms since start of fMRI scan
           % column 10 of datafile is the image
           tmp8 = datafile{8};
           tmp8 = tmp8(2:end);
           % make tmp8 an array of onsets
           onsets = zeros(length(tmp8),1);
           for i = 1:length(tmp8)
               onsets(i) = str2num(tmp8{i});
           end
           tmp10 = datafile{10};
           tmp10 = tmp10(2:end); % this is the list of images
           index_LC     = [];
           index_HC     = [];
           index_NF     = [];
           index_rest   = [];
           index_rating = [];
           for i = 1:length(tmp10)
               if strfind(tmp10{i},'Low') == 1
                  index_LC = [index_LC i];
               elseif strfind(tmp10{i},'High') == 1
                  index_HC = [index_HC i];
               elseif strfind(tmp10{i},'Neutral') == 1
                  index_NF = [index_NF i];
               elseif strfind(tmp10{i},'flow') == 1
                  index_rest = [index_rest i];
               end
           end
           % since 8 images were shown, the index_LC, index_HC and index_NF
           % contain always blocks of 8 subsequent images. We only need the
           % onset of the block of images.
           onsets_blocks{1} = onsets(index_rest)/1000; % in seconds
           onsets_blocks{2} = onsets(index_HC(1:8:end))/1000; % in seconds
           onsets_blocks{3} = onsets(index_LC(1:8:end))/1000; % in seconds
           onsets_blocks{4} = onsets(index_NF(1:8:end))/1000; % in seconds
           tmp1 = onsets(index_HC(8:8:end))/1000 + 3; % onset of the rating block after HC
           tmp1_durations = 12*ones(length(tmp1),1);
           tmp2 = onsets(index_LC(8:8:end))/1000 + 3; % onset of the rating block after LC 
           tmp2_durations = 12*ones(length(tmp2),1);
           tmp3 = onsets(index_NF(8:8:end))/1000 + 3; % onset of the rating block after NF
           tmp3_durations = 6*ones(length(tmp3),1);
           tmp = [tmp1; tmp2; tmp3];
           tmp_durations = [tmp1_durations; tmp2_durations; tmp3_durations];
           [onsets_blocks{5},listorder] = sort(tmp,'ascend');
           durations{1} = 24; % rest period of 24s
           durations{2} = 24; % 8 HC images x 3s
           durations{3} = 24; % 8 LC images x 3s
           durations{4} = 24; % 8 NF images x 3s        
           durations{5} = tmp_durations(listorder);

           for i = 1:nr_conditions12
               matlabbatch{1}.spm.stats.fmri_spec.sess(batch_index).cond(i).name = conditions_name12{i};
               matlabbatch{1}.spm.stats.fmri_spec.sess(batch_index).cond(i).onset = onsets_blocks{i};
               matlabbatch{1}.spm.stats.fmri_spec.sess(batch_index).cond(i).duration = durations{i};
               matlabbatch{1}.spm.stats.fmri_spec.sess(batch_index).cond(i).tmod = 0;
               matlabbatch{1}.spm.stats.fmri_spec.sess(batch_index).cond(i).pmod = struct('name', {}, 'param', {}, 'poly', {});
               matlabbatch{1}.spm.stats.fmri_spec.sess(batch_index).cond(i).orth = 1;
           end
           matlabbatch{1}.spm.stats.fmri_spec.sess(batch_index).multi_reg = {''};
           matlabbatch{1}.spm.stats.fmri_spec.sess(batch_index).hpf = 128;
           matlabbatch{1}.spm.stats.fmri_spec.sess(batch_index).multi = {''};
    
           % add motion regressors
           %----------------------
           cd(datadir)
           filename_rp = fullfile(datadir,['rp_' subjectcode '_func_con_run' num2str(sess) '.txt']);
           rp_data = load(filename_rp);
           matlabbatch{1}.spm.stats.fmri_spec.sess(batch_index).regress(1).name = 'translation X';
           matlabbatch{1}.spm.stats.fmri_spec.sess(batch_index).regress(2).name = 'translation y';
           matlabbatch{1}.spm.stats.fmri_spec.sess(batch_index).regress(3).name = 'translation z';
           matlabbatch{1}.spm.stats.fmri_spec.sess(batch_index).regress(4).name = 'pitch';
           matlabbatch{1}.spm.stats.fmri_spec.sess(batch_index).regress(5).name = 'roll';
           matlabbatch{1}.spm.stats.fmri_spec.sess(batch_index).regress(6).name = 'yaw';
           matlabbatch{1}.spm.stats.fmri_spec.sess(batch_index).regress(1).val = rp_data(1:nr_volumes1,1);
           matlabbatch{1}.spm.stats.fmri_spec.sess(batch_index).regress(2).val = rp_data(1:nr_volumes1,2);
           matlabbatch{1}.spm.stats.fmri_spec.sess(batch_index).regress(3).val = rp_data(1:nr_volumes1,3);
           matlabbatch{1}.spm.stats.fmri_spec.sess(batch_index).regress(4).val = rp_data(1:nr_volumes1,4);
           matlabbatch{1}.spm.stats.fmri_spec.sess(batch_index).regress(5).val = rp_data(1:nr_volumes1,5);
           matlabbatch{1}.spm.stats.fmri_spec.sess(batch_index).regress(6).val = rp_data(1:nr_volumes1,6);
        end
    end
    
    % run 3 and 4 maybe missing in some subjects
    for sess = 3:4
        filename = fullfile(datadir,[template_prefix subjectcode '_func_con_run' num2str(sess) '.nii']);
        if exist(filename,'file') == 2
           batch_index = batch_index+1;
           if sess == 3
              go3 = 1;
           elseif sess == 4
              go4 = 1;
           end
           % read the scans
           filelist = {};
           for i = 1:nr_volumes2
               filelist{i,1} = fullfile(datadir,[template_prefix subjectcode '_func_con_run' num2str(sess) '.nii,' num2str(i)]);
           end
           matlabbatch{1}.spm.stats.fmri_spec.sess(batch_index).scans = filelist;
             
           % read the logfile
           cd(logdir)
           tmp = dir([subjectcode '_con_run' num2str(sess) '*' template_log2]);        
           fid   = fopen(char(tmp(1).name),'r'); 
           tmp          = textscan(fid,'%s',1,'delimiter','\n'); 
           creationtime = textscan(fid,'%s',1,'delimiter','\n'); 
           tmp          = textscan(fid,'%s',1,'delimiter','\n'); 
           headerstring = textscan(fid,'%s',13,'delimiter','\t'); 
           % headerstring should now contain the header of the log file
           %    'Subject'
           %    'Trial'
           %    'Event Type'
           %    'Code'
           %    'Time'
           %    'TTime'
           %    'Uncertainty'
           %    'Duration'
           %    'Uncertainty'
           %    'ReqTime'
           %    'ReqDur'
           %    'Stim Type'
           %    'Pair Index'
           tmp          = textscan(fid,'%s',1,'delimiter','\n');  % empty line        
           % read now the data part
%            datafile = textscan(fid,'%s%d%s%s%f%f%f%f%f%f%s%s%f','delimiter','\t');
           datafile = textscan(fid,'%s%d%s%s%f%*[^\n]','delimiter','\t');
           fclose(fid);

            %    condition          code in logfile
            %    'cue C0'           open_circle
            %    'cue C2'           one_line
            %    'cue C10'          two_lines
            %    'cor after C10'    feedback_10SP AND [(right_arrow followed by response 12) OR (left_arrow followed by response 11)]
            %    'incor after C10'  feedback_10SP AND [(right_arrow followed by response 12) OR (left_arrow followed by response 11)]
            %    'cor after C2'     feedback_2SP AND [(right_arrow followed by response 12) OR (left_arrow followed by response 11)] 
            %    'incor after C2'   feedback_2SP AND [(right_arrow followed by response 11) OR (left_arrow followed by response 12)]
            %    'after C0'         feedback_0SP AND [(right_arrow followed by response 11 or 12) OR (left_arrow followed by response 11 or 12)]
  
           % get the onsets for each condition
           timing = datafile{5}/10000; % in seconds
           code = datafile{4};
           index_cue_C0     = [];
           index_cue_C2     = [];
           index_cue_C10    = [];
           index_feedback_C10 = [];
           index_feedback_C2  = [];
           index_feedback_C0  = [];
           index_Rarrow  = [];
           index_Larrow  = [];
           index_responseL = [];
           index_responseR = [];
           index_scan = []; % all the pulses from the scanner
           for i = 1:length(code)
               if strcmp(code{i},'open_circle')
                  index_cue_C0 = [index_cue_C0 i]; 
               elseif strcmp(code{i},'one_line')
                  index_cue_C2 = [index_cue_C2 i]; 
               elseif strcmp(code{i},'two_lines')
                  index_cue_C10 = [index_cue_C10 i]; 
               elseif strcmp(code{i},'feedback_10SP')
                  index_feedback_C10 = [index_feedback_C10 i];
               elseif strcmp(code{i},'feedback_2SP')
                  index_feedback_C2 = [index_feedback_C2 i];
               elseif strcmp(code{i},'feedback_0SP')
                  index_feedback_C0 = [index_feedback_C0 i];
               elseif strcmp(code{i},'right_arrow')
                  index_Rarrow = [index_Rarrow i];
               elseif strcmp(code{i},'left_arrow')
                  index_Larrow = [index_Larrow i];
               elseif strcmp(code{i},'11')
                  index_responseL = [index_responseL i];
               elseif strcmp(code{i},'12')
                  index_responseR = [index_responseR i];
               elseif strcmp(code{i},'99')
                  index_scan = [index_scan i];
               end
           end
           index_trials_start_cue = sort([index_cue_C0 index_cue_C2 index_cue_C10]);
           % determine the onsets taken into account that the logfile
           % starts when the dummy scan started. So we have to correct the
           % timing with respect to the onset of the fifth scan.
           %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
           clear onsets
           onsets{1} = timing(index_cue_C0)-timing(index_scan(5)); 
           onsets{2} = timing(index_cue_C2)-timing(index_scan(5));
           onsets{3} = timing(index_cue_C10)-timing(index_scan(5));
           timing_feedback_0SP = timing(index_feedback_C0)-timing(index_scan(5));
           timing_feedback_2SP = timing(index_feedback_C2)-timing(index_scan(5));
           timing_feedback_10SP = timing(index_feedback_C10)-timing(index_scan(5));
           timing_Larrow = timing(index_Larrow)-timing(index_scan(5));
           timing_Rarrow = timing(index_Rarrow)-timing(index_scan(5));
           timing_Lresponse = timing(index_responseL)-timing(index_scan(5));
           timing_Rresponse = timing(index_responseR)-timing(index_scan(5));
           timing_feedback = timing(index_trials_start_cue)-timing(index_scan(5))+0.750+3+1; % 0.750s cue display, 3s delay, 1s target

           % determine which one corresponds to a correct or incorrect response
           onsets{4} = []; % 'cor after C10' feedback_10SP AND [(right_arrow followed by response 12) OR (left_arrow followed by response 11)]
           onsets{5} = []; % 'incor after C10' feedback_10SP AND [(right_arrow followed by response 12) OR (left_arrow followed by response 11)]
           onsets{6} = []; % 'cor after C2' feedback_2SP AND [(right_arrow followed by response 12) OR (left_arrow followed by response 11)]
           onsets{7} = []; % 'incor after C2' feedback_2SP AND [(right_arrow followed by response 11) OR (left_arrow followed by response 12)]
           onsets{8} = []; % 'after C0' feedback_0SP AND [(right_arrow followed by response 11 or 12) OR (left_arrow followed by response 11 or 12)]
           for i = 1:length(timing_feedback_10SP)
               time0 = timing_feedback_10SP(i);
               La = intersect(find(timing_Larrow < time0),find(timing_Larrow > time0-2.0)); % actually -1.5s but due to small differences we can miss one. 
               Ra = intersect(find(timing_Rarrow < time0),find(timing_Rarrow > time0-2.0)); % actually -1.5s but due to small differences we can miss one. 
               Lr = intersect(find(timing_Lresponse < time0),find(timing_Lresponse > time0-2.0)); % actually -1.5s but due to small differences we can miss one. 
               Rr = intersect(find(timing_Rresponse < time0),find(timing_Rresponse > time0-2.0)); % actually -1.5s but due to small differences we can miss one. 
               if (length(La) == 1 && length(Lr) == 1) || (length(Ra) == 1 && length(Rr) == 1) % correct response
                  onsets{4} = [onsets{4}; time0];
               end
           end
           for i = 1:length(timing_feedback_2SP)
               time0 = timing_feedback_2SP(i);
               La = intersect(find(timing_Larrow < time0),find(timing_Larrow > time0-2.0)); % actually -1.5s but due to small differences we can miss one. 
               Ra = intersect(find(timing_Rarrow < time0),find(timing_Rarrow > time0-2.0)); % actually -1.5s but due to small differences we can miss one. 
               Lr = intersect(find(timing_Lresponse < time0),find(timing_Lresponse > time0-2.0)); % actually -1.5s but due to small differences we can miss one. 
               Rr = intersect(find(timing_Rresponse < time0),find(timing_Rresponse > time0-2.0)); % actually -1.5s but due to small differences we can miss one. 
               if (length(La) == 1 && length(Lr) == 1) || (length(Ra) == 1 && length(Rr) == 1) % correct response
                  % now check if this is a win trial or a no win trial
                  onsets{6} = [onsets{6}; time0];
               end
           end
           for i = 1:length(timing_feedback_0SP)
               time0 = timing_feedback_0SP(i);
               La = intersect(find(timing_Larrow < time0),find(timing_Larrow > time0-2.0)); % actually -1.5s but due to small differences we can miss one. 
               Ra = intersect(find(timing_Rarrow < time0),find(timing_Rarrow > time0-2.0)); % actually -1.5s but due to small differences we can miss one. 
               Lr = intersect(find(timing_Lresponse < time0),find(timing_Lresponse > time0-2.0)); % actually -1.5s but due to small differences we can miss one. 
               Rr = intersect(find(timing_Rresponse < time0),find(timing_Rresponse > time0-2.0)); % actually -1.5s but due to small differences we can miss one. 
               if (length(La) == 1 && length(Lr) == 1) || (length(Ra) == 1 && length(Rr) == 1) % correct response
                  % now check if this is a win trial or a no win trial
                  index_trial = find(abs(timing_feedback-time0) == min(abs(timing_feedback-time0)));
                  if order_trials(index_trial) == 2
                     onsets{5} = [onsets{5}; time0];
                  elseif order_trials(index_trial) == 4
                     onsets{7} = [onsets{7}; time0];
                  elseif order_trials(index_trial) == 5
                     onsets{8} = [onsets{8}; time0];
                  else
                     fprintf('Something is wrong: inconsistencies found for feedback_C0 \n');
                  end
               elseif (length(La) == 1 && length(Rr) == 1) || (length(Ra) == 1 && length(Lr) == 1)  % incorrect response
                  % now check if this is a win trial or a no win trial
                  index_trial = find(abs(timing_feedback-time0) == min(abs(timing_feedback-time0)));
                  if order_trials(index_trial) == 2
                     onsets{5} = [onsets{5}; time0];
                  elseif order_trials(index_trial) == 4
                     onsets{7} = [onsets{7}; time0];
                  elseif order_trials(index_trial) == 5
                     onsets{8} = [onsets{8}; time0];
                  else
                     fprintf('Something is wrong: inconsistencies found for feedback_C0 \n');
                  end
               end
           end
           
           for i = 1:nr_conditions34
               matlabbatch{1}.spm.stats.fmri_spec.sess(batch_index).cond(i).name = conditions_name34{i};
               matlabbatch{1}.spm.stats.fmri_spec.sess(batch_index).cond(i).onset = onsets{i};
               matlabbatch{1}.spm.stats.fmri_spec.sess(batch_index).cond(i).duration = 0; % all modelled as events
               matlabbatch{1}.spm.stats.fmri_spec.sess(batch_index).cond(i).tmod = 0;
               matlabbatch{1}.spm.stats.fmri_spec.sess(batch_index).cond(i).pmod = struct('name', {}, 'param', {}, 'poly', {});
               matlabbatch{1}.spm.stats.fmri_spec.sess(batch_index).cond(i).orth = 1;
           end
           matlabbatch{1}.spm.stats.fmri_spec.sess(batch_index).multi_reg = {''};
           matlabbatch{1}.spm.stats.fmri_spec.sess(batch_index).hpf = 128;
           matlabbatch{1}.spm.stats.fmri_spec.sess(batch_index).multi = {''};
    
           % add motion regressors
           %----------------------
           cd(datadir)
           filename_rp = fullfile(datadir,['rp_' subjectcode '_func_con_run' num2str(sess) '.txt']);
           rp_data = load(filename_rp);
           matlabbatch{1}.spm.stats.fmri_spec.sess(batch_index).regress(1).name = 'translation X';
           matlabbatch{1}.spm.stats.fmri_spec.sess(batch_index).regress(2).name = 'translation y';
           matlabbatch{1}.spm.stats.fmri_spec.sess(batch_index).regress(3).name = 'translation z';
           matlabbatch{1}.spm.stats.fmri_spec.sess(batch_index).regress(4).name = 'pitch';
           matlabbatch{1}.spm.stats.fmri_spec.sess(batch_index).regress(5).name = 'roll';
           matlabbatch{1}.spm.stats.fmri_spec.sess(batch_index).regress(6).name = 'yaw';
           matlabbatch{1}.spm.stats.fmri_spec.sess(batch_index).regress(1).val = rp_data(1:nr_volumes2,1);
           matlabbatch{1}.spm.stats.fmri_spec.sess(batch_index).regress(2).val = rp_data(1:nr_volumes2,2);
           matlabbatch{1}.spm.stats.fmri_spec.sess(batch_index).regress(3).val = rp_data(1:nr_volumes2,3);
           matlabbatch{1}.spm.stats.fmri_spec.sess(batch_index).regress(4).val = rp_data(1:nr_volumes2,4);
           matlabbatch{1}.spm.stats.fmri_spec.sess(batch_index).regress(5).val = rp_data(1:nr_volumes2,5);
           matlabbatch{1}.spm.stats.fmri_spec.sess(batch_index).regress(6).val = rp_data(1:nr_volumes2,6);
        end
    end
    
    % SECOND BATCH ESTIMATION
    %------------------------
    matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('fMRI model specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
    matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
    matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
    
    % THIRD BATCH CONTRASTS
    %----------------------
    matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
    % contrasts for run1 and run2 (taken into account 6 motion regressors
    % and 5 conditions)
    % 1) HC - NF    0  1  0 -1 0   0 0 0 0 0 0   0  1  0 -1 0   0 0 0 0 0 0
    % 2) NF - HC    0 -1  0  1 0   0 0 0 0 0 0   0 -1  0  1 0   0 0 0 0 0 0
    % 3) LC - NF    0  0  1 -1 0   0 0 0 0 0 0   0  0  1 -1 0   0 0 0 0 0 0
    % 4) NF - LC    0  0 -1  1 0   0 0 0 0 0 0   0  0 -1  1 0   0 0 0 0 0 0
    % 5) HC - LC    0  1 -1  0 0   0 0 0 0 0 0   0  1 -1  0 0   0 0 0 0 0 0
    % 6) LC - HC    0 -1  1  0 0   0 0 0 0 0 0   0 -1  1  0 0   0 0 0 0 0 0
    %
    % contrasts for run3 and run4 ((taken into account 6 motion regressors
    % and 8 conditions, and 11 regressors each for run1 and run2)
    % 7)  cue C10 - cue C0    zeros(1,22) -1  0  1 0 0 0 0 0   0 0 0 0 0 0   -1  0  1 0 0 0 0 0   0 0 0 0 0 0              
    % 8)  cue C0 - cue C10    zeros(1,22)  1  0 -1 0 0 0 0 0   0 0 0 0 0 0    1  0 -1 0 0 0 0 0   0 0 0 0 0 0
    % 9)  cue C2 - cue C0     zeros(1,22) -1  1  0 0 0 0 0 0   0 0 0 0 0 0   -1  1  0 0 0 0 0 0   0 0 0 0 0 0
    % 10) cue C0 - cue C2     zeros(1,22)  1 -1  0 0 0 0 0 0   0 0 0 0 0 0    1 -1  0 0 0 0 0 0   0 0 0 0 0 0
    % 11) cue C10 - cue C2    zeros(1,22)  0 -1  1 0 0 0 0 0   0 0 0 0 0 0    0 -1  1 0 0 0 0 0   0 0 0 0 0 0
    % 12) cue C2 - cue C10    zeros(1,22)  0  1 -1 0 0 0 0 0   0 0 0 0 0 0    0  1 -1 0 0 0 0 0   0 0 0 0 0 0
    % 13) feedback_C10_win - feedback_C0        zeros(1,22) 0 0 0  1 0  0 0 -1   0 0 0 0 0 0    0 0 0  1 0  0 0 -1   0 0 0 0 0 0
    % 14) feedback_C0 - feedback_C10_win        zeros(1,22) 0 0 0 -1 0  0 0  1   0 0 0 0 0 0    0 0 0 -1 0  0 0  1   0 0 0 0 0 0
    % 15) feedback_C2_win - feedback_C0         zeros(1,22) 0 0 0  0 0  1 0 -1   0 0 0 0 0 0    0 0 0  0 0  1 0 -1   0 0 0 0 0 0
    % 16) feedback_C0 - feedback_C2_win         zeros(1,22) 0 0 0  0 0 -1 0  1   0 0 0 0 0 0    0 0 0  0 0 -1 0  1   0 0 0 0 0 0
    % 17) feedback_C10_win - feedback_C2_win    zeros(1,22) 0 0 0  1 0 -1 0  0   0 0 0 0 0 0    0 0 0  1 0 -1 0  0   0 0 0 0 0 0
    % 18) feedback_C2_win - feedback_C10_win    zeros(1,22) 0 0 0 -1 0  1 0  0   0 0 0 0 0 0    0 0 0 -1 0  1 0  0   0 0 0 0 0 0
    if go1 + go2 == 2
       matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'HC - NF';
       matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = [0  1  0 -1 0   0 0 0 0 0 0   0  1  0 -1 0   0 0 0 0 0 0];
       matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
    
       matlabbatch{3}.spm.stats.con.consess{2}.tcon.name = 'NF - HC';
       matlabbatch{3}.spm.stats.con.consess{2}.tcon.weights = [0 -1  0  1 0   0 0 0 0 0 0   0 -1  0  1 0   0 0 0 0 0 0];
       matlabbatch{3}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
    
       matlabbatch{3}.spm.stats.con.consess{3}.tcon.name = 'LC - NF';
       matlabbatch{3}.spm.stats.con.consess{3}.tcon.weights = [0  0  1 -1 0   0 0 0 0 0 0   0  0  1 -1 0   0 0 0 0 0 0];
       matlabbatch{3}.spm.stats.con.consess{3}.tcon.sessrep = 'none';
    
       matlabbatch{3}.spm.stats.con.consess{4}.tcon.name = 'NF - LC';
       matlabbatch{3}.spm.stats.con.consess{4}.tcon.weights = [0  0 -1  1 0   0 0 0 0 0 0   0  0 -1  1 0   0 0 0 0 0 0];
       matlabbatch{3}.spm.stats.con.consess{4}.tcon.sessrep = 'none';
    
       matlabbatch{3}.spm.stats.con.consess{5}.tcon.name = 'HC - LC';
       matlabbatch{3}.spm.stats.con.consess{5}.tcon.weights = [0  1 -1  0 0   0 0 0 0 0 0   0  1 -1  0 0   0 0 0 0 0 0];
       matlabbatch{3}.spm.stats.con.consess{5}.tcon.sessrep = 'none';
    
       matlabbatch{3}.spm.stats.con.consess{6}.tcon.name = 'LC - HC';
       matlabbatch{3}.spm.stats.con.consess{6}.tcon.weights = [0 -1  1  0 0   0 0 0 0 0 0   0 -1  1  0 0   0 0 0 0 0 0];
       matlabbatch{3}.spm.stats.con.consess{6}.tcon.sessrep = 'none';
    elseif go1 + go2 == 1
       matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'HC - NF';
       matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = [0  1  0 -1 0   0 0 0 0 0 0];
       matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
    
       matlabbatch{3}.spm.stats.con.consess{2}.tcon.name = 'NF - HC';
       matlabbatch{3}.spm.stats.con.consess{2}.tcon.weights = [0 -1  0  1 0   0 0 0 0 0 0];
       matlabbatch{3}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
    
       matlabbatch{3}.spm.stats.con.consess{3}.tcon.name = 'LC - NF';
       matlabbatch{3}.spm.stats.con.consess{3}.tcon.weights = [0  0  1 -1 0   0 0 0 0 0 0];
       matlabbatch{3}.spm.stats.con.consess{3}.tcon.sessrep = 'none';
    
       matlabbatch{3}.spm.stats.con.consess{4}.tcon.name = 'NF - LC';
       matlabbatch{3}.spm.stats.con.consess{4}.tcon.weights = [0  0 -1  1 0   0 0 0 0 0 0];
       matlabbatch{3}.spm.stats.con.consess{4}.tcon.sessrep = 'none';
    
       matlabbatch{3}.spm.stats.con.consess{5}.tcon.name = 'HC - LC';
       matlabbatch{3}.spm.stats.con.consess{5}.tcon.weights = [0  1 -1  0 0   0 0 0 0 0 0];
       matlabbatch{3}.spm.stats.con.consess{5}.tcon.sessrep = 'none';
    
       matlabbatch{3}.spm.stats.con.consess{6}.tcon.name = 'LC - HC';
       matlabbatch{3}.spm.stats.con.consess{6}.tcon.weights = [0 -1  1  0 0   0 0 0 0 0 0];
       matlabbatch{3}.spm.stats.con.consess{6}.tcon.sessrep = 'none';
    end        
        
    if go1+go2+go3+go4 == 4
       matlabbatch{3}.spm.stats.con.consess{7}.tcon.name = 'cue C10 - cue C0';
       matlabbatch{3}.spm.stats.con.consess{7}.tcon.weights = [zeros(1,22) -1  0  1 0 0 0 0 0   0 0 0 0 0 0   -1  0  1 0 0 0 0 0   0 0 0 0 0 0];
       matlabbatch{3}.spm.stats.con.consess{7}.tcon.sessrep = 'none';

       matlabbatch{3}.spm.stats.con.consess{8}.tcon.name = 'cue C0 - cue C10';
       matlabbatch{3}.spm.stats.con.consess{8}.tcon.weights = [zeros(1,22)  1  0 -1 0 0 0 0 0   0 0 0 0 0 0    1  0 -1 0 0 0 0 0   0 0 0 0 0 0];
       matlabbatch{3}.spm.stats.con.consess{8}.tcon.sessrep = 'none';

       matlabbatch{3}.spm.stats.con.consess{9}.tcon.name = 'cue C2 - cue C0';
       matlabbatch{3}.spm.stats.con.consess{9}.tcon.weights = [zeros(1,22) -1  1  0 0 0 0 0 0   0 0 0 0 0 0   -1  1  0 0 0 0 0 0   0 0 0 0 0 0];
       matlabbatch{3}.spm.stats.con.consess{9}.tcon.sessrep = 'none';

       matlabbatch{3}.spm.stats.con.consess{10}.tcon.name = 'cue C0 - cue C2';
       matlabbatch{3}.spm.stats.con.consess{10}.tcon.weights = [zeros(1,22)  1 -1  0 0 0 0 0 0   0 0 0 0 0 0    1 -1  0 0 0 0 0 0   0 0 0 0 0 0];
       matlabbatch{3}.spm.stats.con.consess{10}.tcon.sessrep = 'none';

       matlabbatch{3}.spm.stats.con.consess{11}.tcon.name = 'cue C10 - cue C2';
       matlabbatch{3}.spm.stats.con.consess{11}.tcon.weights = [zeros(1,22)  0 -1  1 0 0 0 0 0   0 0 0 0 0 0    0 -1  1 0 0 0 0 0   0 0 0 0 0 0];
       matlabbatch{3}.spm.stats.con.consess{11}.tcon.sessrep = 'none';

       matlabbatch{3}.spm.stats.con.consess{12}.tcon.name = 'cue C2 - cue C10';
       matlabbatch{3}.spm.stats.con.consess{12}.tcon.weights = [zeros(1,22)  0  1 -1 0 0 0 0 0   0 0 0 0 0 0    0  1 -1 0 0 0 0 0   0 0 0 0 0 0];
       matlabbatch{3}.spm.stats.con.consess{12}.tcon.sessrep = 'none';
    
       matlabbatch{3}.spm.stats.con.consess{13}.tcon.name = 'feedback_C10_win - feedback_C0';
       matlabbatch{3}.spm.stats.con.consess{13}.tcon.weights = [zeros(1,22) 0 0 0  1 0  0 0 -1   0 0 0 0 0 0    0 0 0  1 0  0 0 -1   0 0 0 0 0 0];
       matlabbatch{3}.spm.stats.con.consess{13}.tcon.sessrep = 'none';

       matlabbatch{3}.spm.stats.con.consess{14}.tcon.name = 'feedback_C0 - feedback_C10_win';
       matlabbatch{3}.spm.stats.con.consess{14}.tcon.weights = [zeros(1,22) 0 0 0 -1 0  0 0  1   0 0 0 0 0 0    0 0 0 -1 0  0 0  1   0 0 0 0 0 0];
       matlabbatch{3}.spm.stats.con.consess{14}.tcon.sessrep = 'none';

       matlabbatch{3}.spm.stats.con.consess{15}.tcon.name = 'feedback_C2_win - feedback_C0';
       matlabbatch{3}.spm.stats.con.consess{15}.tcon.weights = [zeros(1,22) 0 0 0  0 0  1 0 -1   0 0 0 0 0 0    0 0 0  0 0  1 0 -1   0 0 0 0 0 0];
       matlabbatch{3}.spm.stats.con.consess{15}.tcon.sessrep = 'none';
    
       matlabbatch{3}.spm.stats.con.consess{16}.tcon.name = 'feedback_C0 - feedback_C2_win';
       matlabbatch{3}.spm.stats.con.consess{16}.tcon.weights = [zeros(1,22) 0 0 0  0 0 -1 0  1   0 0 0 0 0 0    0 0 0  0 0 -1 0  1   0 0 0 0 0 0];
       matlabbatch{3}.spm.stats.con.consess{16}.tcon.sessrep = 'none';
    
       matlabbatch{3}.spm.stats.con.consess{17}.tcon.name = 'feedback_C10_win - feedback_C2_win';
       matlabbatch{3}.spm.stats.con.consess{17}.tcon.weights = [zeros(1,22) 0 0 0  1 0 -1 0  0   0 0 0 0 0 0    0 0 0  1 0 -1 0  0   0 0 0 0 0 0];
       matlabbatch{3}.spm.stats.con.consess{17}.tcon.sessrep = 'none';
    
       matlabbatch{3}.spm.stats.con.consess{18}.tcon.name = 'feedback_C2_win - feedback_C10_win';
       matlabbatch{3}.spm.stats.con.consess{18}.tcon.weights = [zeros(1,22) 0 0 0 -1 0  1 0  0   0 0 0 0 0 0    0 0 0 -1 0  1 0  0   0 0 0 0 0 0];
       matlabbatch{3}.spm.stats.con.consess{18}.tcon.sessrep = 'none';
    elseif go1+go2 == 2 && go3+go4 == 1
       matlabbatch{3}.spm.stats.con.consess{7}.tcon.name = 'cue C10 - cue C0';
       matlabbatch{3}.spm.stats.con.consess{7}.tcon.weights = [zeros(1,22) -1  0  1 0 0 0 0 0   0 0 0 0 0 0];
       matlabbatch{3}.spm.stats.con.consess{7}.tcon.sessrep = 'none';

       matlabbatch{3}.spm.stats.con.consess{8}.tcon.name = 'cue C0 - cue C10';
       matlabbatch{3}.spm.stats.con.consess{8}.tcon.weights = [zeros(1,22)  1  0 -1 0 0 0 0 0   0 0 0 0 0 0];
       matlabbatch{3}.spm.stats.con.consess{8}.tcon.sessrep = 'none';

       matlabbatch{3}.spm.stats.con.consess{9}.tcon.name = 'cue C2 - cue C0';
       matlabbatch{3}.spm.stats.con.consess{9}.tcon.weights = [zeros(1,22) -1  1  0 0 0 0 0 0   0 0 0 0 0 0];
       matlabbatch{3}.spm.stats.con.consess{9}.tcon.sessrep = 'none';

       matlabbatch{3}.spm.stats.con.consess{10}.tcon.name = 'cue C0 - cue C2';
       matlabbatch{3}.spm.stats.con.consess{10}.tcon.weights = [zeros(1,22)  1 -1  0 0 0 0 0 0   0 0 0 0 0 0];
       matlabbatch{3}.spm.stats.con.consess{10}.tcon.sessrep = 'none';

       matlabbatch{3}.spm.stats.con.consess{11}.tcon.name = 'cue C10 - cue C2';
       matlabbatch{3}.spm.stats.con.consess{11}.tcon.weights = [zeros(1,22)  0 -1  1 0 0 0 0 0   0 0 0 0 0 0];
       matlabbatch{3}.spm.stats.con.consess{11}.tcon.sessrep = 'none';

       matlabbatch{3}.spm.stats.con.consess{12}.tcon.name = 'cue C2 - cue C10';
       matlabbatch{3}.spm.stats.con.consess{12}.tcon.weights = [zeros(1,22)  0  1 -1 0 0 0 0 0   0 0 0 0 0 0];
       matlabbatch{3}.spm.stats.con.consess{12}.tcon.sessrep = 'none';
    
       matlabbatch{3}.spm.stats.con.consess{13}.tcon.name = 'feedback_C10_win - feedback_C0';
       matlabbatch{3}.spm.stats.con.consess{13}.tcon.weights = [zeros(1,22) 0 0 0  1 0  0 0 -1   0 0 0 0 0 0];
       matlabbatch{3}.spm.stats.con.consess{13}.tcon.sessrep = 'none';

       matlabbatch{3}.spm.stats.con.consess{14}.tcon.name = 'feedback_C0 - feedback_C10_win';
       matlabbatch{3}.spm.stats.con.consess{14}.tcon.weights = [zeros(1,22) 0 0 0 -1 0  0 0  1   0 0 0 0 0 0];
       matlabbatch{3}.spm.stats.con.consess{14}.tcon.sessrep = 'none';

       matlabbatch{3}.spm.stats.con.consess{15}.tcon.name = 'feedback_C2_win - feedback_C0';
       matlabbatch{3}.spm.stats.con.consess{15}.tcon.weights = [zeros(1,22) 0 0 0  0 0  1 0 -1   0 0 0 0 0 0];
       matlabbatch{3}.spm.stats.con.consess{15}.tcon.sessrep = 'none';
    
       matlabbatch{3}.spm.stats.con.consess{16}.tcon.name = 'feedback_C0 - feedback_C2_win';
       matlabbatch{3}.spm.stats.con.consess{16}.tcon.weights = [zeros(1,22) 0 0 0  0 0 -1 0  1   0 0 0 0 0 0];
       matlabbatch{3}.spm.stats.con.consess{16}.tcon.sessrep = 'none';
    
       matlabbatch{3}.spm.stats.con.consess{17}.tcon.name = 'feedback_C10_win - feedback_C2_win';
       matlabbatch{3}.spm.stats.con.consess{17}.tcon.weights = [zeros(1,22) 0 0 0  1 0 -1 0  0   0 0 0 0 0 0];
       matlabbatch{3}.spm.stats.con.consess{17}.tcon.sessrep = 'none';
    
       matlabbatch{3}.spm.stats.con.consess{18}.tcon.name = 'feedback_C2_win - feedback_C10_win';
       matlabbatch{3}.spm.stats.con.consess{18}.tcon.weights = [zeros(1,22) 0 0 0 -1 0  1 0  0   0 0 0 0 0 0];
       matlabbatch{3}.spm.stats.con.consess{18}.tcon.sessrep = 'none';
    elseif go1+go2 == 1 && go3+go4 == 2
       matlabbatch{3}.spm.stats.con.consess{7}.tcon.name = 'cue C10 - cue C0';
       matlabbatch{3}.spm.stats.con.consess{7}.tcon.weights = [zeros(1,11) -1  0  1 0 0 0 0 0   0 0 0 0 0 0   -1  0  1 0 0 0 0 0   0 0 0 0 0 0];
       matlabbatch{3}.spm.stats.con.consess{7}.tcon.sessrep = 'none';

       matlabbatch{3}.spm.stats.con.consess{8}.tcon.name = 'cue C0 - cue C10';
       matlabbatch{3}.spm.stats.con.consess{8}.tcon.weights = [zeros(1,11)  1  0 -1 0 0 0 0 0   0 0 0 0 0 0    1  0 -1 0 0 0 0 0   0 0 0 0 0 0];
       matlabbatch{3}.spm.stats.con.consess{8}.tcon.sessrep = 'none';

       matlabbatch{3}.spm.stats.con.consess{9}.tcon.name = 'cue C2 - cue C0';
       matlabbatch{3}.spm.stats.con.consess{9}.tcon.weights = [zeros(1,11) -1  1  0 0 0 0 0 0   0 0 0 0 0 0   -1  1  0 0 0 0 0 0   0 0 0 0 0 0];
       matlabbatch{3}.spm.stats.con.consess{9}.tcon.sessrep = 'none';

       matlabbatch{3}.spm.stats.con.consess{10}.tcon.name = 'cue C0 - cue C2';
       matlabbatch{3}.spm.stats.con.consess{10}.tcon.weights = [zeros(1,11)  1 -1  0 0 0 0 0 0   0 0 0 0 0 0    1 -1  0 0 0 0 0 0   0 0 0 0 0 0];
       matlabbatch{3}.spm.stats.con.consess{10}.tcon.sessrep = 'none';

       matlabbatch{3}.spm.stats.con.consess{11}.tcon.name = 'cue C10 - cue C2';
       matlabbatch{3}.spm.stats.con.consess{11}.tcon.weights = [zeros(1,11)  0 -1  1 0 0 0 0 0   0 0 0 0 0 0    0 -1  1 0 0 0 0 0   0 0 0 0 0 0];
       matlabbatch{3}.spm.stats.con.consess{11}.tcon.sessrep = 'none';

       matlabbatch{3}.spm.stats.con.consess{12}.tcon.name = 'cue C2 - cue C10';
       matlabbatch{3}.spm.stats.con.consess{12}.tcon.weights = [zeros(1,11)  0  1 -1 0 0 0 0 0   0 0 0 0 0 0    0  1 -1 0 0 0 0 0   0 0 0 0 0 0];
       matlabbatch{3}.spm.stats.con.consess{12}.tcon.sessrep = 'none';
    
       matlabbatch{3}.spm.stats.con.consess{13}.tcon.name = 'feedback_C10_win - feedback_C0';
       matlabbatch{3}.spm.stats.con.consess{13}.tcon.weights = [zeros(1,11) 0 0 0  1 0  0 0 -1   0 0 0 0 0 0    0 0 0  1 0  0 0 -1   0 0 0 0 0 0];
       matlabbatch{3}.spm.stats.con.consess{13}.tcon.sessrep = 'none';

       matlabbatch{3}.spm.stats.con.consess{14}.tcon.name = 'feedback_C0 - feedback_C10_win';
       matlabbatch{3}.spm.stats.con.consess{14}.tcon.weights = [zeros(1,11) 0 0 0 -1 0  0 0  1   0 0 0 0 0 0    0 0 0 -1 0  0 0  1   0 0 0 0 0 0];
       matlabbatch{3}.spm.stats.con.consess{14}.tcon.sessrep = 'none';

       matlabbatch{3}.spm.stats.con.consess{15}.tcon.name = 'feedback_C2_win - feedback_C0';
       matlabbatch{3}.spm.stats.con.consess{15}.tcon.weights = [zeros(1,11) 0 0 0  0 0  1 0 -1   0 0 0 0 0 0    0 0 0  0 0  1 0 -1   0 0 0 0 0 0];
       matlabbatch{3}.spm.stats.con.consess{15}.tcon.sessrep = 'none';
    
       matlabbatch{3}.spm.stats.con.consess{16}.tcon.name = 'feedback_C0 - feedback_C2_win';
       matlabbatch{3}.spm.stats.con.consess{16}.tcon.weights = [zeros(1,11) 0 0 0  0 0 -1 0  1   0 0 0 0 0 0    0 0 0  0 0 -1 0  1   0 0 0 0 0 0];
       matlabbatch{3}.spm.stats.con.consess{16}.tcon.sessrep = 'none';
    
       matlabbatch{3}.spm.stats.con.consess{17}.tcon.name = 'feedback_C10_win - feedback_C2_win';
       matlabbatch{3}.spm.stats.con.consess{17}.tcon.weights = [zeros(1,11) 0 0 0  1 0 -1 0  0   0 0 0 0 0 0    0 0 0  1 0 -1 0  0   0 0 0 0 0 0];
       matlabbatch{3}.spm.stats.con.consess{17}.tcon.sessrep = 'none';
    
       matlabbatch{3}.spm.stats.con.consess{18}.tcon.name = 'feedback_C2_win - feedback_C10_win';
       matlabbatch{3}.spm.stats.con.consess{18}.tcon.weights = [zeros(1,11) 0 0 0 -1 0  1 0  0   0 0 0 0 0 0    0 0 0 -1 0  1 0  0   0 0 0 0 0 0];
       matlabbatch{3}.spm.stats.con.consess{18}.tcon.sessrep = 'none';
    elseif go1+go2 == 1 && go3+go4 == 1
       matlabbatch{3}.spm.stats.con.consess{7}.tcon.name = 'cue C10 - cue C0';
       matlabbatch{3}.spm.stats.con.consess{7}.tcon.weights = [zeros(1,11) -1  0  1 0 0 0 0 0   0 0 0 0 0 0];
       matlabbatch{3}.spm.stats.con.consess{7}.tcon.sessrep = 'none';

       matlabbatch{3}.spm.stats.con.consess{8}.tcon.name = 'cue C0 - cue C10';
       matlabbatch{3}.spm.stats.con.consess{8}.tcon.weights = [zeros(1,11)  1  0 -1 0 0 0 0 0   0 0 0 0 0 0];
       matlabbatch{3}.spm.stats.con.consess{8}.tcon.sessrep = 'none';

       matlabbatch{3}.spm.stats.con.consess{9}.tcon.name = 'cue C2 - cue C0';
       matlabbatch{3}.spm.stats.con.consess{9}.tcon.weights = [zeros(1,11) -1  1  0 0 0 0 0 0   0 0 0 0 0 0];
       matlabbatch{3}.spm.stats.con.consess{9}.tcon.sessrep = 'none';

       matlabbatch{3}.spm.stats.con.consess{10}.tcon.name = 'cue C0 - cue C2';
       matlabbatch{3}.spm.stats.con.consess{10}.tcon.weights = [zeros(1,11)  1 -1  0 0 0 0 0 0   0 0 0 0 0 0];
       matlabbatch{3}.spm.stats.con.consess{10}.tcon.sessrep = 'none';

       matlabbatch{3}.spm.stats.con.consess{11}.tcon.name = 'cue C10 - cue C2';
       matlabbatch{3}.spm.stats.con.consess{11}.tcon.weights = [zeros(1,11)  0 -1  1 0 0 0 0 0   0 0 0 0 0 0];
       matlabbatch{3}.spm.stats.con.consess{11}.tcon.sessrep = 'none';

       matlabbatch{3}.spm.stats.con.consess{12}.tcon.name = 'cue C2 - cue C10';
       matlabbatch{3}.spm.stats.con.consess{12}.tcon.weights = [zeros(1,11)  0  1 -1 0 0 0 0 0   0 0 0 0 0 0];
       matlabbatch{3}.spm.stats.con.consess{12}.tcon.sessrep = 'none';
    
       matlabbatch{3}.spm.stats.con.consess{13}.tcon.name = 'feedback_C10_win - feedback_C0';
       matlabbatch{3}.spm.stats.con.consess{13}.tcon.weights = [zeros(1,11) 0 0 0  1 0  0 0 -1   0 0 0 0 0 0];
       matlabbatch{3}.spm.stats.con.consess{13}.tcon.sessrep = 'none';

       matlabbatch{3}.spm.stats.con.consess{14}.tcon.name = 'feedback_C0 - feedback_C10_win';
       matlabbatch{3}.spm.stats.con.consess{14}.tcon.weights = [zeros(1,11) 0 0 0 -1 0  0 0  1   0 0 0 0 0 0];
       matlabbatch{3}.spm.stats.con.consess{14}.tcon.sessrep = 'none';

       matlabbatch{3}.spm.stats.con.consess{15}.tcon.name = 'feedback_C2_win - feedback_C0';
       matlabbatch{3}.spm.stats.con.consess{15}.tcon.weights = [zeros(1,11) 0 0 0  0 0  1 0 -1   0 0 0 0 0 0];
       matlabbatch{3}.spm.stats.con.consess{15}.tcon.sessrep = 'none';
    
       matlabbatch{3}.spm.stats.con.consess{16}.tcon.name = 'feedback_C0 - feedback_C2_win';
       matlabbatch{3}.spm.stats.con.consess{16}.tcon.weights = [zeros(1,11) 0 0 0  0 0 -1 0  1   0 0 0 0 0 0];
       matlabbatch{3}.spm.stats.con.consess{16}.tcon.sessrep = 'none';
    
       matlabbatch{3}.spm.stats.con.consess{17}.tcon.name = 'feedback_C10_win - feedback_C2_win';
       matlabbatch{3}.spm.stats.con.consess{17}.tcon.weights = [zeros(1,11) 0 0 0  1 0 -1 0  0   0 0 0 0 0 0];
       matlabbatch{3}.spm.stats.con.consess{17}.tcon.sessrep = 'none';
    
       matlabbatch{3}.spm.stats.con.consess{18}.tcon.name = 'feedback_C2_win - feedback_C10_win';
       matlabbatch{3}.spm.stats.con.consess{18}.tcon.weights = [zeros(1,11) 0 0 0 -1 0  1 0  0   0 0 0 0 0 0];
       matlabbatch{3}.spm.stats.con.consess{18}.tcon.sessrep = 'none';
    end 
    matlabbatch{3}.spm.stats.con.delete = 0;
   
    
    cd(outputdir_con)
    save batch_stats matlabbatch;  
    spm_jobman('initcfg')
    spm_jobman('run',matlabbatch);%run batch
    
    % FOR THE BIT STUDY
    %++++++++++++++++++
        
    go1 = 0; % for run1
    go2 = 0; % for run2
    go3 = 0; % for run3
    go4 = 0; % for run4
    batch_index = 0;

    % FIRST BATCH MODEL SPECIFICATION
    %--------------------------------
    % create empty structures
    matlabbatch = struct([]);
    
    matlabbatch{1}.spm.stats.fmri_spec.dir = {outputdir_bit};
    matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs';
    matlabbatch{1}.spm.stats.fmri_spec.timing.RT = TR;
    matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 16;
    matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 1;
    for sess = 1:2
        filename = fullfile(datadir,[template_prefix subjectcode '_func_bit_run' num2str(sess) '.nii']);
        if exist(filename,'file') == 2
           batch_index = batch_index+1;
           if sess == 1
              go1 = 1;
           elseif sess == 2
              go2 = 1;
           end
       % read the scans
           filelist = {};
           for i = 1:nr_volumes1
               filelist{i,1} = fullfile(datadir,[template_prefix subjectcode '_func_bit_run' num2str(sess) '.nii,' num2str(i)]);
           end
           matlabbatch{1}.spm.stats.fmri_spec.sess(batch_index).scans = filelist;
             
           % read the logfile
           cd(logdir)
           tmp = dir([subjectcode '_bit_run' num2str(sess) '*' template_log1]);        
           fid   = fopen(char(tmp(1).name),'r'); 
           datafile = textscan(fid,'%s %s %s %s %s %s %s %s %s %s %*[^\n]','delimiter','\t');
           fclose(fid);
           % first row of datafile is a header and should be neglected
           % column 8 of datafile is the start time in ms since start of fMRI scan
           % column 9 of datafile is the end time in ms since start of fMRI scan
           % column 10 of datafile is the image
           tmp8 = datafile{8};
           tmp8 = tmp8(2:end);
           % make tmp8 an array of onsets
           onsets = zeros(length(tmp8),1);
           for i = 1:length(tmp8)
               onsets(i) = str2num(tmp8{i});
           end
           tmp10 = datafile{10};
           tmp10 = tmp10(2:end); % this is the list of images
           index_LC     = [];
           index_HC     = [];
           index_NF     = [];
           index_rest   = [];
           index_rating = [];
           for i = 1:length(tmp10)
               if strfind(tmp10{i},'Low') == 1
                  index_LC = [index_LC i];
               elseif strfind(tmp10{i},'High') == 1
                  index_HC = [index_HC i];
               elseif strfind(tmp10{i},'Neutral') == 1
                  index_NF = [index_NF i];
               elseif strfind(tmp10{i},'flow') == 1
                  index_rest = [index_rest i];
               end
           end
           % since 8 images were shown, the index_LC, index_HC and index_NF
           % bittain always blocks of 8 subsequent images. We only need the
           % onset of the block of images.
           onsets_blocks{1} = onsets(index_rest)/1000; % in seconds
           onsets_blocks{2} = onsets(index_HC(1:8:end))/1000; % in seconds
           onsets_blocks{3} = onsets(index_LC(1:8:end))/1000; % in seconds
           onsets_blocks{4} = onsets(index_NF(1:8:end))/1000; % in seconds
           tmp1 = onsets(index_HC(8:8:end))/1000 + 3; % onset of the rating block after HC
           tmp1_durations = 12*ones(length(tmp1),1);
           tmp2 = onsets(index_LC(8:8:end))/1000 + 3; % onset of the rating block after LC 
           tmp2_durations = 12*ones(length(tmp2),1);
           tmp3 = onsets(index_NF(8:8:end))/1000 + 3; % onset of the rating block after NF
           tmp3_durations = 6*ones(length(tmp3),1);
           tmp = [tmp1; tmp2; tmp3];
           tmp_durations = [tmp1_durations; tmp2_durations; tmp3_durations];
           [onsets_blocks{5},listorder] = sort(tmp,'ascend');
           durations{1} = 24; % rest period of 24s
           durations{2} = 24; % 8 HC images x 3s
           durations{3} = 24; % 8 LC images x 3s
           durations{4} = 24; % 8 NF images x 3s        
           durations{5} = tmp_durations(listorder);

           for i = 1:nr_conditions12
               matlabbatch{1}.spm.stats.fmri_spec.sess(batch_index).cond(i).name = conditions_name12{i};
               matlabbatch{1}.spm.stats.fmri_spec.sess(batch_index).cond(i).onset = onsets_blocks{i};
               matlabbatch{1}.spm.stats.fmri_spec.sess(batch_index).cond(i).duration = durations{i};
               matlabbatch{1}.spm.stats.fmri_spec.sess(batch_index).cond(i).tmod = 0;
               matlabbatch{1}.spm.stats.fmri_spec.sess(batch_index).cond(i).pmod = struct('name', {}, 'param', {}, 'poly', {});
               matlabbatch{1}.spm.stats.fmri_spec.sess(batch_index).cond(i).orth = 1;
           end
           matlabbatch{1}.spm.stats.fmri_spec.sess(batch_index).multi_reg = {''};
           matlabbatch{1}.spm.stats.fmri_spec.sess(batch_index).hpf = 128;
           matlabbatch{1}.spm.stats.fmri_spec.sess(batch_index).multi = {''};
    
           % add motion regressors
           %----------------------
           cd(datadir)
           filename_rp = fullfile(datadir,['rp_' subjectcode '_func_bit_run' num2str(sess) '.txt']);
           rp_data = load(filename_rp);
           matlabbatch{1}.spm.stats.fmri_spec.sess(batch_index).regress(1).name = 'translation X';
           matlabbatch{1}.spm.stats.fmri_spec.sess(batch_index).regress(2).name = 'translation y';
           matlabbatch{1}.spm.stats.fmri_spec.sess(batch_index).regress(3).name = 'translation z';
           matlabbatch{1}.spm.stats.fmri_spec.sess(batch_index).regress(4).name = 'pitch';
           matlabbatch{1}.spm.stats.fmri_spec.sess(batch_index).regress(5).name = 'roll';
           matlabbatch{1}.spm.stats.fmri_spec.sess(batch_index).regress(6).name = 'yaw';
           matlabbatch{1}.spm.stats.fmri_spec.sess(batch_index).regress(1).val = rp_data(1:nr_volumes1,1);
           matlabbatch{1}.spm.stats.fmri_spec.sess(batch_index).regress(2).val = rp_data(1:nr_volumes1,2);
           matlabbatch{1}.spm.stats.fmri_spec.sess(batch_index).regress(3).val = rp_data(1:nr_volumes1,3);
           matlabbatch{1}.spm.stats.fmri_spec.sess(batch_index).regress(4).val = rp_data(1:nr_volumes1,4);
           matlabbatch{1}.spm.stats.fmri_spec.sess(batch_index).regress(5).val = rp_data(1:nr_volumes1,5);
           matlabbatch{1}.spm.stats.fmri_spec.sess(batch_index).regress(6).val = rp_data(1:nr_volumes1,6);
        end
    end
    
    % run 3 and 4 maybe missing in some subjects
    for sess = 3:4
        filename = fullfile(datadir,[template_prefix subjectcode '_func_bit_run' num2str(sess) '.nii']);
        if exist(filename,'file') == 2
           batch_index = batch_index+1;
           if sess == 3
              go3 = 1;
           elseif sess == 4
              go4 = 1;
           end
           % read the scans
           filelist = {};
           for i = 1:nr_volumes2
               filelist{i,1} = fullfile(datadir,[template_prefix subjectcode '_func_bit_run' num2str(sess) '.nii,' num2str(i)]);
           end
           matlabbatch{1}.spm.stats.fmri_spec.sess(batch_index).scans = filelist;
             
           % read the logfile
           cd(logdir)
           tmp = dir([subjectcode '_bit_run' num2str(sess) '*' template_log2]);        
           fid   = fopen(char(tmp(1).name),'r'); 
           tmp          = textscan(fid,'%s',1,'delimiter','\n'); 
           creationtime = textscan(fid,'%s',1,'delimiter','\n'); 
           tmp          = textscan(fid,'%s',1,'delimiter','\n'); 
           headerstring = textscan(fid,'%s',13,'delimiter','\t'); 
           % headerstring should now contain the header of the log file
           %    'Subject'
           %    'Trial'
           %    'Event Type'
           %    'Code'
           %    'Time'
           %    'TTime'
           %    'Uncertainty'
           %    'Duration'
           %    'Uncertainty'
           %    'ReqTime'
           %    'ReqDur'
           %    'Stim Type'
           %    'Pair Index'
           tmp          = textscan(fid,'%s',1,'delimiter','\n');  % empty line        
           % read now the data part
%            datafile = textscan(fid,'%s%d%s%s%f%f%f%f%f%f%s%s%f','delimiter','\t');
           datafile = textscan(fid,'%s%d%s%s%f%*[^\n]','delimiter','\t');
           fclose(fid);

            %    condition          code in logfile
            %    'cue C0'           open_circle
            %    'cue C2'           one_line
            %    'cue C10'          two_lines
            %    'cor after C10'    feedback_10SP AND [(right_arrow followed by response 12) OR (left_arrow followed by response 11)]
            %    'incor after C10'  feedback_10SP AND [(right_arrow followed by response 12) OR (left_arrow followed by response 11)]
            %    'cor after C2'     feedback_2SP AND [(right_arrow followed by response 12) OR (left_arrow followed by response 11)] 
            %    'incor after C2'   feedback_2SP AND [(right_arrow followed by response 11) OR (left_arrow followed by response 12)]
            %    'after C0'         feedback_0SP AND [(right_arrow followed by response 11 or 12) OR (left_arrow followed by response 11 or 12)]
  
           % get the onsets for each condition
           timing = datafile{5}/10000; % in seconds
           code = datafile{4};
           index_cue_C0     = [];
           index_cue_C2     = [];
           index_cue_C10    = [];
           index_feedback_C10 = [];
           index_feedback_C2  = [];
           index_feedback_C0  = [];
           index_Rarrow  = [];
           index_Larrow  = [];
           index_responseL = [];
           index_responseR = [];
           index_scan = []; % all the pulses from the scanner
           for i = 1:length(code)
               if strcmp(code{i},'open_circle')
                  index_cue_C0 = [index_cue_C0 i]; 
               elseif strcmp(code{i},'one_line')
                  index_cue_C2 = [index_cue_C2 i]; 
               elseif strcmp(code{i},'two_lines')
                  index_cue_C10 = [index_cue_C10 i]; 
               elseif strcmp(code{i},'feedback_10SP')
                  index_feedback_C10 = [index_feedback_C10 i];
               elseif strcmp(code{i},'feedback_2SP')
                  index_feedback_C2 = [index_feedback_C2 i];
               elseif strcmp(code{i},'feedback_0SP')
                  index_feedback_C0 = [index_feedback_C0 i];
               elseif strcmp(code{i},'right_arrow')
                  index_Rarrow = [index_Rarrow i];
               elseif strcmp(code{i},'left_arrow')
                  index_Larrow = [index_Larrow i];
               elseif strcmp(code{i},'11')
                  index_responseL = [index_responseL i];
               elseif strcmp(code{i},'12')
                  index_responseR = [index_responseR i];
               elseif strcmp(code{i},'99')
                  index_scan = [index_scan i];
               end
           end
           index_trials_start_cue = sort([index_cue_C0 index_cue_C2 index_cue_C10]);
           % determine the onsets taken into account that the logfile
           % starts when the dummy scan started. So we have to correct the
           % timing with respect to the onset of the fifth scan.
           %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
           clear onsets
           onsets{1} = timing(index_cue_C0)-timing(index_scan(5)); 
           onsets{2} = timing(index_cue_C2)-timing(index_scan(5));
           onsets{3} = timing(index_cue_C10)-timing(index_scan(5));
           timing_feedback_0SP = timing(index_feedback_C0)-timing(index_scan(5));
           timing_feedback_2SP = timing(index_feedback_C2)-timing(index_scan(5));
           timing_feedback_10SP = timing(index_feedback_C10)-timing(index_scan(5));
           timing_Larrow = timing(index_Larrow)-timing(index_scan(5));
           timing_Rarrow = timing(index_Rarrow)-timing(index_scan(5));
           timing_Lresponse = timing(index_responseL)-timing(index_scan(5));
           timing_Rresponse = timing(index_responseR)-timing(index_scan(5));
           timing_feedback = timing(index_trials_start_cue)-timing(index_scan(5))+0.750+3+1; % 0.750s cue display, 3s delay, 1s target

           % determine which one corresponds to a correct or incorrect response
           onsets{4} = []; % 'cor after C10' feedback_10SP AND [(right_arrow followed by response 12) OR (left_arrow followed by response 11)]
           onsets{5} = []; % 'incor after C10' feedback_10SP AND [(right_arrow followed by response 12) OR (left_arrow followed by response 11)]
           onsets{6} = []; % 'cor after C2' feedback_2SP AND [(right_arrow followed by response 12) OR (left_arrow followed by response 11)]
           onsets{7} = []; % 'incor after C2' feedback_2SP AND [(right_arrow followed by response 11) OR (left_arrow followed by response 12)]
           onsets{8} = []; % 'after C0' feedback_0SP AND [(right_arrow followed by response 11 or 12) OR (left_arrow followed by response 11 or 12)]
           for i = 1:length(timing_feedback_10SP)
               time0 = timing_feedback_10SP(i);
               La = intersect(find(timing_Larrow < time0),find(timing_Larrow > time0-2.0)); % actually -1.5s but due to small differences we can miss one. 
               Ra = intersect(find(timing_Rarrow < time0),find(timing_Rarrow > time0-2.0)); % actually -1.5s but due to small differences we can miss one. 
               Lr = intersect(find(timing_Lresponse < time0),find(timing_Lresponse > time0-2.0)); % actually -1.5s but due to small differences we can miss one. 
               Rr = intersect(find(timing_Rresponse < time0),find(timing_Rresponse > time0-2.0)); % actually -1.5s but due to small differences we can miss one. 
               if (length(La) == 1 && length(Lr) == 1) || (length(Ra) == 1 && length(Rr) == 1) % correct response
                  onsets{4} = [onsets{4}; time0];
               end
           end
           for i = 1:length(timing_feedback_2SP)
               time0 = timing_feedback_2SP(i);
               La = intersect(find(timing_Larrow < time0),find(timing_Larrow > time0-2.0)); % actually -1.5s but due to small differences we can miss one. 
               Ra = intersect(find(timing_Rarrow < time0),find(timing_Rarrow > time0-2.0)); % actually -1.5s but due to small differences we can miss one. 
               Lr = intersect(find(timing_Lresponse < time0),find(timing_Lresponse > time0-2.0)); % actually -1.5s but due to small differences we can miss one. 
               Rr = intersect(find(timing_Rresponse < time0),find(timing_Rresponse > time0-2.0)); % actually -1.5s but due to small differences we can miss one. 
               if (length(La) == 1 && length(Lr) == 1) || (length(Ra) == 1 && length(Rr) == 1) % correct response
                  % now check if this is a win trial or a no win trial
                  onsets{6} = [onsets{6}; time0];
               end
           end
           for i = 1:length(timing_feedback_0SP)
               time0 = timing_feedback_0SP(i);
               La = intersect(find(timing_Larrow < time0),find(timing_Larrow > time0-2.0)); % actually -1.5s but due to small differences we can miss one. 
               Ra = intersect(find(timing_Rarrow < time0),find(timing_Rarrow > time0-2.0)); % actually -1.5s but due to small differences we can miss one. 
               Lr = intersect(find(timing_Lresponse < time0),find(timing_Lresponse > time0-2.0)); % actually -1.5s but due to small differences we can miss one. 
               Rr = intersect(find(timing_Rresponse < time0),find(timing_Rresponse > time0-2.0)); % actually -1.5s but due to small differences we can miss one. 
               if (length(La) == 1 && length(Lr) == 1) || (length(Ra) == 1 && length(Rr) == 1) % correct response
                  % now check if this is a win trial or a no win trial
                  index_trial = find(abs(timing_feedback-time0) == min(abs(timing_feedback-time0)));
                  if order_trials(index_trial) == 2
                     onsets{5} = [onsets{5}; time0];
                  elseif order_trials(index_trial) == 4
                     onsets{7} = [onsets{7}; time0];
                  elseif order_trials(index_trial) == 5
                     onsets{8} = [onsets{8}; time0];
                  else
                     fprintf('Something is wrong: inconsistencies found for feedback_C0 \n');
                  end
               elseif (length(La) == 1 && length(Rr) == 1) || (length(Ra) == 1 && length(Lr) == 1)  % incorrect response
                  % now check if this is a win trial or a no win trial
                  index_trial = find(abs(timing_feedback-time0) == min(abs(timing_feedback-time0)));
                  if order_trials(index_trial) == 2
                     onsets{5} = [onsets{5}; time0];
                  elseif order_trials(index_trial) == 4
                     onsets{7} = [onsets{7}; time0];
                  elseif order_trials(index_trial) == 5
                     onsets{8} = [onsets{8}; time0];
                  else
                     fprintf('Something is wrong: inconsistencies found for feedback_C0 \n');
                  end
               end
           end
           
           for i = 1:nr_conditions34
               matlabbatch{1}.spm.stats.fmri_spec.sess(batch_index).cond(i).name = conditions_name34{i};
               matlabbatch{1}.spm.stats.fmri_spec.sess(batch_index).cond(i).onset = onsets{i};
               matlabbatch{1}.spm.stats.fmri_spec.sess(batch_index).cond(i).duration = 0; % all modelled as events
               matlabbatch{1}.spm.stats.fmri_spec.sess(batch_index).cond(i).tmod = 0;
               matlabbatch{1}.spm.stats.fmri_spec.sess(batch_index).cond(i).pmod = struct('name', {}, 'param', {}, 'poly', {});
               matlabbatch{1}.spm.stats.fmri_spec.sess(batch_index).cond(i).orth = 1;
           end
           matlabbatch{1}.spm.stats.fmri_spec.sess(batch_index).multi_reg = {''};
           matlabbatch{1}.spm.stats.fmri_spec.sess(batch_index).hpf = 128;
           matlabbatch{1}.spm.stats.fmri_spec.sess(batch_index).multi = {''};
    
           % add motion regressors
           %----------------------
           cd(datadir)
           filename_rp = fullfile(datadir,['rp_' subjectcode '_func_bit_run' num2str(sess) '.txt']);
           rp_data = load(filename_rp);
           matlabbatch{1}.spm.stats.fmri_spec.sess(batch_index).regress(1).name = 'translation X';
           matlabbatch{1}.spm.stats.fmri_spec.sess(batch_index).regress(2).name = 'translation y';
           matlabbatch{1}.spm.stats.fmri_spec.sess(batch_index).regress(3).name = 'translation z';
           matlabbatch{1}.spm.stats.fmri_spec.sess(batch_index).regress(4).name = 'pitch';
           matlabbatch{1}.spm.stats.fmri_spec.sess(batch_index).regress(5).name = 'roll';
           matlabbatch{1}.spm.stats.fmri_spec.sess(batch_index).regress(6).name = 'yaw';
           matlabbatch{1}.spm.stats.fmri_spec.sess(batch_index).regress(1).val = rp_data(1:nr_volumes2,1);
           matlabbatch{1}.spm.stats.fmri_spec.sess(batch_index).regress(2).val = rp_data(1:nr_volumes2,2);
           matlabbatch{1}.spm.stats.fmri_spec.sess(batch_index).regress(3).val = rp_data(1:nr_volumes2,3);
           matlabbatch{1}.spm.stats.fmri_spec.sess(batch_index).regress(4).val = rp_data(1:nr_volumes2,4);
           matlabbatch{1}.spm.stats.fmri_spec.sess(batch_index).regress(5).val = rp_data(1:nr_volumes2,5);
           matlabbatch{1}.spm.stats.fmri_spec.sess(batch_index).regress(6).val = rp_data(1:nr_volumes2,6);
        end
    end
    
    % SECOND BATCH ESTIMATION
    %------------------------
    matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('fMRI model specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
    matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
    matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
    
    % THIRD BATCH CONTRASTS
    %----------------------
    matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
    % contrasts for run1 and run2 (taken into account 6 motion regressors
    % and 5 conditions)
    % 1) HC - NF    0  1  0 -1 0   0 0 0 0 0 0   0  1  0 -1 0   0 0 0 0 0 0
    % 2) NF - HC    0 -1  0  1 0   0 0 0 0 0 0   0 -1  0  1 0   0 0 0 0 0 0
    % 3) LC - NF    0  0  1 -1 0   0 0 0 0 0 0   0  0  1 -1 0   0 0 0 0 0 0
    % 4) NF - LC    0  0 -1  1 0   0 0 0 0 0 0   0  0 -1  1 0   0 0 0 0 0 0
    % 5) HC - LC    0  1 -1  0 0   0 0 0 0 0 0   0  1 -1  0 0   0 0 0 0 0 0
    % 6) LC - HC    0 -1  1  0 0   0 0 0 0 0 0   0 -1  1  0 0   0 0 0 0 0 0
    %
    % contrasts for run3 and run4 ((taken into account 6 motion regressors
    % and 8 conditions, and 11 regressors each for run1 and run2)
    % 7)  cue C10 - cue C0    zeros(1,22) -1  0  1 0 0 0 0 0   0 0 0 0 0 0   -1  0  1 0 0 0 0 0   0 0 0 0 0 0              
    % 8)  cue C0 - cue C10    zeros(1,22)  1  0 -1 0 0 0 0 0   0 0 0 0 0 0    1  0 -1 0 0 0 0 0   0 0 0 0 0 0
    % 9)  cue C2 - cue C0     zeros(1,22) -1  1  0 0 0 0 0 0   0 0 0 0 0 0   -1  1  0 0 0 0 0 0   0 0 0 0 0 0
    % 10) cue C0 - cue C2     zeros(1,22)  1 -1  0 0 0 0 0 0   0 0 0 0 0 0    1 -1  0 0 0 0 0 0   0 0 0 0 0 0
    % 11) cue C10 - cue C2    zeros(1,22)  0 -1  1 0 0 0 0 0   0 0 0 0 0 0    0 -1  1 0 0 0 0 0   0 0 0 0 0 0
    % 12) cue C2 - cue C10    zeros(1,22)  0  1 -1 0 0 0 0 0   0 0 0 0 0 0    0  1 -1 0 0 0 0 0   0 0 0 0 0 0
    % 13) feedback_C10_win - feedback_C0        zeros(1,22) 0 0 0  1 0  0 0 -1   0 0 0 0 0 0    0 0 0  1 0  0 0 -1   0 0 0 0 0 0
    % 14) feedback_C0 - feedback_C10_win        zeros(1,22) 0 0 0 -1 0  0 0  1   0 0 0 0 0 0    0 0 0 -1 0  0 0  1   0 0 0 0 0 0
    % 15) feedback_C2_win - feedback_C0         zeros(1,22) 0 0 0  0 0  1 0 -1   0 0 0 0 0 0    0 0 0  0 0  1 0 -1   0 0 0 0 0 0
    % 16) feedback_C0 - feedback_C2_win         zeros(1,22) 0 0 0  0 0 -1 0  1   0 0 0 0 0 0    0 0 0  0 0 -1 0  1   0 0 0 0 0 0
    % 17) feedback_C10_win - feedback_C2_win    zeros(1,22) 0 0 0  1 0 -1 0  0   0 0 0 0 0 0    0 0 0  1 0 -1 0  0   0 0 0 0 0 0
    % 18) feedback_C2_win - feedback_C10_win    zeros(1,22) 0 0 0 -1 0  1 0  0   0 0 0 0 0 0    0 0 0 -1 0  1 0  0   0 0 0 0 0 0
    if go1 + go2 == 2
       matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'HC - NF';
       matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = [0  1  0 -1 0   0 0 0 0 0 0   0  1  0 -1 0   0 0 0 0 0 0];
       matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
    
       matlabbatch{3}.spm.stats.con.consess{2}.tcon.name = 'NF - HC';
       matlabbatch{3}.spm.stats.con.consess{2}.tcon.weights = [0 -1  0  1 0   0 0 0 0 0 0   0 -1  0  1 0   0 0 0 0 0 0];
       matlabbatch{3}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
    
       matlabbatch{3}.spm.stats.con.consess{3}.tcon.name = 'LC - NF';
       matlabbatch{3}.spm.stats.con.consess{3}.tcon.weights = [0  0  1 -1 0   0 0 0 0 0 0   0  0  1 -1 0   0 0 0 0 0 0];
       matlabbatch{3}.spm.stats.con.consess{3}.tcon.sessrep = 'none';
    
       matlabbatch{3}.spm.stats.con.consess{4}.tcon.name = 'NF - LC';
       matlabbatch{3}.spm.stats.con.consess{4}.tcon.weights = [0  0 -1  1 0   0 0 0 0 0 0   0  0 -1  1 0   0 0 0 0 0 0];
       matlabbatch{3}.spm.stats.con.consess{4}.tcon.sessrep = 'none';
    
       matlabbatch{3}.spm.stats.con.consess{5}.tcon.name = 'HC - LC';
       matlabbatch{3}.spm.stats.con.consess{5}.tcon.weights = [0  1 -1  0 0   0 0 0 0 0 0   0  1 -1  0 0   0 0 0 0 0 0];
       matlabbatch{3}.spm.stats.con.consess{5}.tcon.sessrep = 'none';
    
       matlabbatch{3}.spm.stats.con.consess{6}.tcon.name = 'LC - HC';
       matlabbatch{3}.spm.stats.con.consess{6}.tcon.weights = [0 -1  1  0 0   0 0 0 0 0 0   0 -1  1  0 0   0 0 0 0 0 0];
       matlabbatch{3}.spm.stats.con.consess{6}.tcon.sessrep = 'none';
    elseif go1 + go2 == 1
       matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'HC - NF';
       matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = [0  1  0 -1 0   0 0 0 0 0 0];
       matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
    
       matlabbatch{3}.spm.stats.con.consess{2}.tcon.name = 'NF - HC';
       matlabbatch{3}.spm.stats.con.consess{2}.tcon.weights = [0 -1  0  1 0   0 0 0 0 0 0];
       matlabbatch{3}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
    
       matlabbatch{3}.spm.stats.con.consess{3}.tcon.name = 'LC - NF';
       matlabbatch{3}.spm.stats.con.consess{3}.tcon.weights = [0  0  1 -1 0   0 0 0 0 0 0];
       matlabbatch{3}.spm.stats.con.consess{3}.tcon.sessrep = 'none';
    
       matlabbatch{3}.spm.stats.con.consess{4}.tcon.name = 'NF - LC';
       matlabbatch{3}.spm.stats.con.consess{4}.tcon.weights = [0  0 -1  1 0   0 0 0 0 0 0];
       matlabbatch{3}.spm.stats.con.consess{4}.tcon.sessrep = 'none';
    
       matlabbatch{3}.spm.stats.con.consess{5}.tcon.name = 'HC - LC';
       matlabbatch{3}.spm.stats.con.consess{5}.tcon.weights = [0  1 -1  0 0   0 0 0 0 0 0];
       matlabbatch{3}.spm.stats.con.consess{5}.tcon.sessrep = 'none';
    
       matlabbatch{3}.spm.stats.con.consess{6}.tcon.name = 'LC - HC';
       matlabbatch{3}.spm.stats.con.consess{6}.tcon.weights = [0 -1  1  0 0   0 0 0 0 0 0];
       matlabbatch{3}.spm.stats.con.consess{6}.tcon.sessrep = 'none';
    end        
        
    if go1+go2+go3+go4 == 4
       matlabbatch{3}.spm.stats.con.consess{7}.tcon.name = 'cue C10 - cue C0';
       matlabbatch{3}.spm.stats.con.consess{7}.tcon.weights = [zeros(1,22) -1  0  1 0 0 0 0 0   0 0 0 0 0 0   -1  0  1 0 0 0 0 0   0 0 0 0 0 0];
       matlabbatch{3}.spm.stats.con.consess{7}.tcon.sessrep = 'none';

       matlabbatch{3}.spm.stats.con.consess{8}.tcon.name = 'cue C0 - cue C10';
       matlabbatch{3}.spm.stats.con.consess{8}.tcon.weights = [zeros(1,22)  1  0 -1 0 0 0 0 0   0 0 0 0 0 0    1  0 -1 0 0 0 0 0   0 0 0 0 0 0];
       matlabbatch{3}.spm.stats.con.consess{8}.tcon.sessrep = 'none';

       matlabbatch{3}.spm.stats.con.consess{9}.tcon.name = 'cue C2 - cue C0';
       matlabbatch{3}.spm.stats.con.consess{9}.tcon.weights = [zeros(1,22) -1  1  0 0 0 0 0 0   0 0 0 0 0 0   -1  1  0 0 0 0 0 0   0 0 0 0 0 0];
       matlabbatch{3}.spm.stats.con.consess{9}.tcon.sessrep = 'none';

       matlabbatch{3}.spm.stats.con.consess{10}.tcon.name = 'cue C0 - cue C2';
       matlabbatch{3}.spm.stats.con.consess{10}.tcon.weights = [zeros(1,22)  1 -1  0 0 0 0 0 0   0 0 0 0 0 0    1 -1  0 0 0 0 0 0   0 0 0 0 0 0];
       matlabbatch{3}.spm.stats.con.consess{10}.tcon.sessrep = 'none';

       matlabbatch{3}.spm.stats.con.consess{11}.tcon.name = 'cue C10 - cue C2';
       matlabbatch{3}.spm.stats.con.consess{11}.tcon.weights = [zeros(1,22)  0 -1  1 0 0 0 0 0   0 0 0 0 0 0    0 -1  1 0 0 0 0 0   0 0 0 0 0 0];
       matlabbatch{3}.spm.stats.con.consess{11}.tcon.sessrep = 'none';

       matlabbatch{3}.spm.stats.con.consess{12}.tcon.name = 'cue C2 - cue C10';
       matlabbatch{3}.spm.stats.con.consess{12}.tcon.weights = [zeros(1,22)  0  1 -1 0 0 0 0 0   0 0 0 0 0 0    0  1 -1 0 0 0 0 0   0 0 0 0 0 0];
       matlabbatch{3}.spm.stats.con.consess{12}.tcon.sessrep = 'none';
    
       matlabbatch{3}.spm.stats.con.consess{13}.tcon.name = 'feedback_C10_win - feedback_C0';
       matlabbatch{3}.spm.stats.con.consess{13}.tcon.weights = [zeros(1,22) 0 0 0  1 0  0 0 -1   0 0 0 0 0 0    0 0 0  1 0  0 0 -1   0 0 0 0 0 0];
       matlabbatch{3}.spm.stats.con.consess{13}.tcon.sessrep = 'none';

       matlabbatch{3}.spm.stats.con.consess{14}.tcon.name = 'feedback_C0 - feedback_C10_win';
       matlabbatch{3}.spm.stats.con.consess{14}.tcon.weights = [zeros(1,22) 0 0 0 -1 0  0 0  1   0 0 0 0 0 0    0 0 0 -1 0  0 0  1   0 0 0 0 0 0];
       matlabbatch{3}.spm.stats.con.consess{14}.tcon.sessrep = 'none';

       matlabbatch{3}.spm.stats.con.consess{15}.tcon.name = 'feedback_C2_win - feedback_C0';
       matlabbatch{3}.spm.stats.con.consess{15}.tcon.weights = [zeros(1,22) 0 0 0  0 0  1 0 -1   0 0 0 0 0 0    0 0 0  0 0  1 0 -1   0 0 0 0 0 0];
       matlabbatch{3}.spm.stats.con.consess{15}.tcon.sessrep = 'none';
    
       matlabbatch{3}.spm.stats.con.consess{16}.tcon.name = 'feedback_C0 - feedback_C2_win';
       matlabbatch{3}.spm.stats.con.consess{16}.tcon.weights = [zeros(1,22) 0 0 0  0 0 -1 0  1   0 0 0 0 0 0    0 0 0  0 0 -1 0  1   0 0 0 0 0 0];
       matlabbatch{3}.spm.stats.con.consess{16}.tcon.sessrep = 'none';
    
       matlabbatch{3}.spm.stats.con.consess{17}.tcon.name = 'feedback_C10_win - feedback_C2_win';
       matlabbatch{3}.spm.stats.con.consess{17}.tcon.weights = [zeros(1,22) 0 0 0  1 0 -1 0  0   0 0 0 0 0 0    0 0 0  1 0 -1 0  0   0 0 0 0 0 0];
       matlabbatch{3}.spm.stats.con.consess{17}.tcon.sessrep = 'none';
    
       matlabbatch{3}.spm.stats.con.consess{18}.tcon.name = 'feedback_C2_win - feedback_C10_win';
       matlabbatch{3}.spm.stats.con.consess{18}.tcon.weights = [zeros(1,22) 0 0 0 -1 0  1 0  0   0 0 0 0 0 0    0 0 0 -1 0  1 0  0   0 0 0 0 0 0];
       matlabbatch{3}.spm.stats.con.consess{18}.tcon.sessrep = 'none';
    elseif go1+go2 == 2 && go3+go4 == 1
       matlabbatch{3}.spm.stats.con.consess{7}.tcon.name = 'cue C10 - cue C0';
       matlabbatch{3}.spm.stats.con.consess{7}.tcon.weights = [zeros(1,22) -1  0  1 0 0 0 0 0   0 0 0 0 0 0];
       matlabbatch{3}.spm.stats.con.consess{7}.tcon.sessrep = 'none';

       matlabbatch{3}.spm.stats.con.consess{8}.tcon.name = 'cue C0 - cue C10';
       matlabbatch{3}.spm.stats.con.consess{8}.tcon.weights = [zeros(1,22)  1  0 -1 0 0 0 0 0   0 0 0 0 0 0];
       matlabbatch{3}.spm.stats.con.consess{8}.tcon.sessrep = 'none';

       matlabbatch{3}.spm.stats.con.consess{9}.tcon.name = 'cue C2 - cue C0';
       matlabbatch{3}.spm.stats.con.consess{9}.tcon.weights = [zeros(1,22) -1  1  0 0 0 0 0 0   0 0 0 0 0 0];
       matlabbatch{3}.spm.stats.con.consess{9}.tcon.sessrep = 'none';

       matlabbatch{3}.spm.stats.con.consess{10}.tcon.name = 'cue C0 - cue C2';
       matlabbatch{3}.spm.stats.con.consess{10}.tcon.weights = [zeros(1,22)  1 -1  0 0 0 0 0 0   0 0 0 0 0 0];
       matlabbatch{3}.spm.stats.con.consess{10}.tcon.sessrep = 'none';

       matlabbatch{3}.spm.stats.con.consess{11}.tcon.name = 'cue C10 - cue C2';
       matlabbatch{3}.spm.stats.con.consess{11}.tcon.weights = [zeros(1,22)  0 -1  1 0 0 0 0 0   0 0 0 0 0 0];
       matlabbatch{3}.spm.stats.con.consess{11}.tcon.sessrep = 'none';

       matlabbatch{3}.spm.stats.con.consess{12}.tcon.name = 'cue C2 - cue C10';
       matlabbatch{3}.spm.stats.con.consess{12}.tcon.weights = [zeros(1,22)  0  1 -1 0 0 0 0 0   0 0 0 0 0 0];
       matlabbatch{3}.spm.stats.con.consess{12}.tcon.sessrep = 'none';
    
       matlabbatch{3}.spm.stats.con.consess{13}.tcon.name = 'feedback_C10_win - feedback_C0';
       matlabbatch{3}.spm.stats.con.consess{13}.tcon.weights = [zeros(1,22) 0 0 0  1 0  0 0 -1   0 0 0 0 0 0];
       matlabbatch{3}.spm.stats.con.consess{13}.tcon.sessrep = 'none';

       matlabbatch{3}.spm.stats.con.consess{14}.tcon.name = 'feedback_C0 - feedback_C10_win';
       matlabbatch{3}.spm.stats.con.consess{14}.tcon.weights = [zeros(1,22) 0 0 0 -1 0  0 0  1   0 0 0 0 0 0];
       matlabbatch{3}.spm.stats.con.consess{14}.tcon.sessrep = 'none';

       matlabbatch{3}.spm.stats.con.consess{15}.tcon.name = 'feedback_C2_win - feedback_C0';
       matlabbatch{3}.spm.stats.con.consess{15}.tcon.weights = [zeros(1,22) 0 0 0  0 0  1 0 -1   0 0 0 0 0 0];
       matlabbatch{3}.spm.stats.con.consess{15}.tcon.sessrep = 'none';
    
       matlabbatch{3}.spm.stats.con.consess{16}.tcon.name = 'feedback_C0 - feedback_C2_win';
       matlabbatch{3}.spm.stats.con.consess{16}.tcon.weights = [zeros(1,22) 0 0 0  0 0 -1 0  1   0 0 0 0 0 0];
       matlabbatch{3}.spm.stats.con.consess{16}.tcon.sessrep = 'none';
    
       matlabbatch{3}.spm.stats.con.consess{17}.tcon.name = 'feedback_C10_win - feedback_C2_win';
       matlabbatch{3}.spm.stats.con.consess{17}.tcon.weights = [zeros(1,22) 0 0 0  1 0 -1 0  0   0 0 0 0 0 0];
       matlabbatch{3}.spm.stats.con.consess{17}.tcon.sessrep = 'none';
    
       matlabbatch{3}.spm.stats.con.consess{18}.tcon.name = 'feedback_C2_win - feedback_C10_win';
       matlabbatch{3}.spm.stats.con.consess{18}.tcon.weights = [zeros(1,22) 0 0 0 -1 0  1 0  0   0 0 0 0 0 0];
       matlabbatch{3}.spm.stats.con.consess{18}.tcon.sessrep = 'none';
    elseif go1+go2 == 1 && go3+go4 == 2
       matlabbatch{3}.spm.stats.con.consess{7}.tcon.name = 'cue C10 - cue C0';
       matlabbatch{3}.spm.stats.con.consess{7}.tcon.weights = [zeros(1,11) -1  0  1 0 0 0 0 0   0 0 0 0 0 0   -1  0  1 0 0 0 0 0   0 0 0 0 0 0];
       matlabbatch{3}.spm.stats.con.consess{7}.tcon.sessrep = 'none';

       matlabbatch{3}.spm.stats.con.consess{8}.tcon.name = 'cue C0 - cue C10';
       matlabbatch{3}.spm.stats.con.consess{8}.tcon.weights = [zeros(1,11)  1  0 -1 0 0 0 0 0   0 0 0 0 0 0    1  0 -1 0 0 0 0 0   0 0 0 0 0 0];
       matlabbatch{3}.spm.stats.con.consess{8}.tcon.sessrep = 'none';

       matlabbatch{3}.spm.stats.con.consess{9}.tcon.name = 'cue C2 - cue C0';
       matlabbatch{3}.spm.stats.con.consess{9}.tcon.weights = [zeros(1,11) -1  1  0 0 0 0 0 0   0 0 0 0 0 0   -1  1  0 0 0 0 0 0   0 0 0 0 0 0];
       matlabbatch{3}.spm.stats.con.consess{9}.tcon.sessrep = 'none';

       matlabbatch{3}.spm.stats.con.consess{10}.tcon.name = 'cue C0 - cue C2';
       matlabbatch{3}.spm.stats.con.consess{10}.tcon.weights = [zeros(1,11)  1 -1  0 0 0 0 0 0   0 0 0 0 0 0    1 -1  0 0 0 0 0 0   0 0 0 0 0 0];
       matlabbatch{3}.spm.stats.con.consess{10}.tcon.sessrep = 'none';

       matlabbatch{3}.spm.stats.con.consess{11}.tcon.name = 'cue C10 - cue C2';
       matlabbatch{3}.spm.stats.con.consess{11}.tcon.weights = [zeros(1,11)  0 -1  1 0 0 0 0 0   0 0 0 0 0 0    0 -1  1 0 0 0 0 0   0 0 0 0 0 0];
       matlabbatch{3}.spm.stats.con.consess{11}.tcon.sessrep = 'none';

       matlabbatch{3}.spm.stats.con.consess{12}.tcon.name = 'cue C2 - cue C10';
       matlabbatch{3}.spm.stats.con.consess{12}.tcon.weights = [zeros(1,11)  0  1 -1 0 0 0 0 0   0 0 0 0 0 0    0  1 -1 0 0 0 0 0   0 0 0 0 0 0];
       matlabbatch{3}.spm.stats.con.consess{12}.tcon.sessrep = 'none';
    
       matlabbatch{3}.spm.stats.con.consess{13}.tcon.name = 'feedback_C10_win - feedback_C0';
       matlabbatch{3}.spm.stats.con.consess{13}.tcon.weights = [zeros(1,11) 0 0 0  1 0  0 0 -1   0 0 0 0 0 0    0 0 0  1 0  0 0 -1   0 0 0 0 0 0];
       matlabbatch{3}.spm.stats.con.consess{13}.tcon.sessrep = 'none';

       matlabbatch{3}.spm.stats.con.consess{14}.tcon.name = 'feedback_C0 - feedback_C10_win';
       matlabbatch{3}.spm.stats.con.consess{14}.tcon.weights = [zeros(1,11) 0 0 0 -1 0  0 0  1   0 0 0 0 0 0    0 0 0 -1 0  0 0  1   0 0 0 0 0 0];
       matlabbatch{3}.spm.stats.con.consess{14}.tcon.sessrep = 'none';

       matlabbatch{3}.spm.stats.con.consess{15}.tcon.name = 'feedback_C2_win - feedback_C0';
       matlabbatch{3}.spm.stats.con.consess{15}.tcon.weights = [zeros(1,11) 0 0 0  0 0  1 0 -1   0 0 0 0 0 0    0 0 0  0 0  1 0 -1   0 0 0 0 0 0];
       matlabbatch{3}.spm.stats.con.consess{15}.tcon.sessrep = 'none';
    
       matlabbatch{3}.spm.stats.con.consess{16}.tcon.name = 'feedback_C0 - feedback_C2_win';
       matlabbatch{3}.spm.stats.con.consess{16}.tcon.weights = [zeros(1,11) 0 0 0  0 0 -1 0  1   0 0 0 0 0 0    0 0 0  0 0 -1 0  1   0 0 0 0 0 0];
       matlabbatch{3}.spm.stats.con.consess{16}.tcon.sessrep = 'none';
    
       matlabbatch{3}.spm.stats.con.consess{17}.tcon.name = 'feedback_C10_win - feedback_C2_win';
       matlabbatch{3}.spm.stats.con.consess{17}.tcon.weights = [zeros(1,11) 0 0 0  1 0 -1 0  0   0 0 0 0 0 0    0 0 0  1 0 -1 0  0   0 0 0 0 0 0];
       matlabbatch{3}.spm.stats.con.consess{17}.tcon.sessrep = 'none';
    
       matlabbatch{3}.spm.stats.con.consess{18}.tcon.name = 'feedback_C2_win - feedback_C10_win';
       matlabbatch{3}.spm.stats.con.consess{18}.tcon.weights = [zeros(1,11) 0 0 0 -1 0  1 0  0   0 0 0 0 0 0    0 0 0 -1 0  1 0  0   0 0 0 0 0 0];
       matlabbatch{3}.spm.stats.con.consess{18}.tcon.sessrep = 'none';
    elseif go1+go2 == 1 && go3+go4 == 1
       matlabbatch{3}.spm.stats.con.consess{7}.tcon.name = 'cue C10 - cue C0';
       matlabbatch{3}.spm.stats.con.consess{7}.tcon.weights = [zeros(1,11) -1  0  1 0 0 0 0 0   0 0 0 0 0 0];
       matlabbatch{3}.spm.stats.con.consess{7}.tcon.sessrep = 'none';

       matlabbatch{3}.spm.stats.con.consess{8}.tcon.name = 'cue C0 - cue C10';
       matlabbatch{3}.spm.stats.con.consess{8}.tcon.weights = [zeros(1,11)  1  0 -1 0 0 0 0 0   0 0 0 0 0 0];
       matlabbatch{3}.spm.stats.con.consess{8}.tcon.sessrep = 'none';

       matlabbatch{3}.spm.stats.con.consess{9}.tcon.name = 'cue C2 - cue C0';
       matlabbatch{3}.spm.stats.con.consess{9}.tcon.weights = [zeros(1,11) -1  1  0 0 0 0 0 0   0 0 0 0 0 0];
       matlabbatch{3}.spm.stats.con.consess{9}.tcon.sessrep = 'none';

       matlabbatch{3}.spm.stats.con.consess{10}.tcon.name = 'cue C0 - cue C2';
       matlabbatch{3}.spm.stats.con.consess{10}.tcon.weights = [zeros(1,11)  1 -1  0 0 0 0 0 0   0 0 0 0 0 0];
       matlabbatch{3}.spm.stats.con.consess{10}.tcon.sessrep = 'none';

       matlabbatch{3}.spm.stats.con.consess{11}.tcon.name = 'cue C10 - cue C2';
       matlabbatch{3}.spm.stats.con.consess{11}.tcon.weights = [zeros(1,11)  0 -1  1 0 0 0 0 0   0 0 0 0 0 0];
       matlabbatch{3}.spm.stats.con.consess{11}.tcon.sessrep = 'none';

       matlabbatch{3}.spm.stats.con.consess{12}.tcon.name = 'cue C2 - cue C10';
       matlabbatch{3}.spm.stats.con.consess{12}.tcon.weights = [zeros(1,11)  0  1 -1 0 0 0 0 0   0 0 0 0 0 0];
       matlabbatch{3}.spm.stats.con.consess{12}.tcon.sessrep = 'none';
    
       matlabbatch{3}.spm.stats.con.consess{13}.tcon.name = 'feedback_C10_win - feedback_C0';
       matlabbatch{3}.spm.stats.con.consess{13}.tcon.weights = [zeros(1,11) 0 0 0  1 0  0 0 -1   0 0 0 0 0 0];
       matlabbatch{3}.spm.stats.con.consess{13}.tcon.sessrep = 'none';

       matlabbatch{3}.spm.stats.con.consess{14}.tcon.name = 'feedback_C0 - feedback_C10_win';
       matlabbatch{3}.spm.stats.con.consess{14}.tcon.weights = [zeros(1,11) 0 0 0 -1 0  0 0  1   0 0 0 0 0 0];
       matlabbatch{3}.spm.stats.con.consess{14}.tcon.sessrep = 'none';

       matlabbatch{3}.spm.stats.con.consess{15}.tcon.name = 'feedback_C2_win - feedback_C0';
       matlabbatch{3}.spm.stats.con.consess{15}.tcon.weights = [zeros(1,11) 0 0 0  0 0  1 0 -1   0 0 0 0 0 0];
       matlabbatch{3}.spm.stats.con.consess{15}.tcon.sessrep = 'none';
    
       matlabbatch{3}.spm.stats.con.consess{16}.tcon.name = 'feedback_C0 - feedback_C2_win';
       matlabbatch{3}.spm.stats.con.consess{16}.tcon.weights = [zeros(1,11) 0 0 0  0 0 -1 0  1   0 0 0 0 0 0];
       matlabbatch{3}.spm.stats.con.consess{16}.tcon.sessrep = 'none';
    
       matlabbatch{3}.spm.stats.con.consess{17}.tcon.name = 'feedback_C10_win - feedback_C2_win';
       matlabbatch{3}.spm.stats.con.consess{17}.tcon.weights = [zeros(1,11) 0 0 0  1 0 -1 0  0   0 0 0 0 0 0];
       matlabbatch{3}.spm.stats.con.consess{17}.tcon.sessrep = 'none';
    
       matlabbatch{3}.spm.stats.con.consess{18}.tcon.name = 'feedback_C2_win - feedback_C10_win';
       matlabbatch{3}.spm.stats.con.consess{18}.tcon.weights = [zeros(1,11) 0 0 0 -1 0  1 0  0   0 0 0 0 0 0];
       matlabbatch{3}.spm.stats.con.consess{18}.tcon.sessrep = 'none';
    end 
    matlabbatch{3}.spm.stats.con.delete = 0;
   
    
    cd(outputdir_bit)
    save batch_stats matlabbatch;  
    spm_jobman('initcfg')
    spm_jobman('run',matlabbatch);%run batch
end