%% bit_rew_prep_s1_write_events_tsv
%
% This script reads logfiles, extracts the onsets, durations, and ratings for
% different conditions, and writes events.tsv files to the BIDS dir for
% each subject
% It also contains an option to write a single phenotype file with
% trial-by-trial ratings for all subjects
% 
% USAGE
%
% Script should be run from the root directory of the superdataset, e.g.
% /data/proj_discoverie
% The script is highly study-specific, as logfiles will vary with design,
% stimulus presentation software used, etc
% Hence, it is provided in LaBGAScore as an example and needs to be
% downloaded and adapted to the code subdataset for your study/project
% This example is from LaBGAS proj_erythritol_4a
% (https://gin.g-node.org/labgas/proj_erythritol_4a)
%
%
% DEPENDENCIES
%
% LaBGAScore Github repo on Matlab path, with subfolders
% https://github.com/labgas/LaBGAScore
%
%
% INPUTS
%
% Presentation .log files in sourcedata dir for each subject
%
%
% OUTPUTS
%
% events.tsv files for each run in BIDS dir for each subject
% phenotype.tsv file in BIDS/phenotype dir (optional)
%
%__________________________________________________________________________
%
% author: Lukas Van Oudenhove
% date:   December, 2021
%
%__________________________________________________________________________
% @(#)% LaBGAScore_prep_s1_write_events_tsv.m         v1.2        
% last modified: 2022/04/21
%
% adapted from LCN12_JULIE_first_level_analysis script by Patrick Dupont
%


%% DEFINE DIRECTORIES, SUBJECTS, RUNS, CONDITIONS, AND IMPORT OPTIONS
%--------------------------------------------------------------------------

bit_rew_prep_s0_define_directories; % lukasvo edited from original LaBGAScore script to enable standalone functioning of proj_ery_4a dataset

subjs2write = {'sub-001'}; % enter subjects separated by comma if you only want to write files for selected subjects e.g. {'sub-01','sub-02'}
pheno_tsv = true; % turn to false if you do not wish to generate a phenotype.tsv file with trial-by-trial ratings; will only work if subjs2write is empty (i.e. when you loop over all your subjects)
pheno_name = 'ratings_online.tsv';

runnames = {'run-1','run-2','run-3','run-4'};
logfilenames = {'*_run1_photo_results.txt','*_run2_photo_results.txt','*_run3-FID_1_final*.log','*_run4-FID_1_final*.log'};
taskname1 = 'food_images_';
taskname2 = 'FID_';
image_labels = {'High*';'Low*';'Neutral*';'flow*'}; % labels for sweet substance delivery in Code var of logfile
% swallow_rinse_labels = {'sucrose_swallowing';'erythritol_swallowing';'sucralose_swallowing';'control_swallowing'}; % labels for swallowing cue presentation after sweet substance delivery in Code var of logfile
% rating_labels = {'Sucrose','Erythritol','Sucralose','Control'}; % labels for start of rating period in Code var of logfile
% fixation_labels = {'fixation_cross','Sucrose fixation cross','Erythritol fixation cross','Sucralose fixation cross','Control fixation cross'}; % labels for fixation cross in Code var of logfile
image_events_interest = {'high calorie','low calorie','neutral'}; % names of events of interest to be written to events.tsv
image_events_nuisance = {'rest','rating'}; % names of nuisance events to be written to events.tsv
% 
varNames = {'version','timestamp','block number','trial number','lost msec','free VRAM','Conditie','time_start','time_end','trial contents'}; % varnames of logfile
selectedVarNames = [8:10]; % varnames we want to use in the script
varTypes = {'char','char','char','char','char','char','char','double','double','char'}; % matlab vartypes to be used when importing log file as table
delimiter = '\t';
dataStartLine = 2; % line on which actual data starts in logfile
extraColRule = 'ignore';

opts = delimitedTextImportOptions('VariableNames',varNames,...
                                'SelectedVariableNames',selectedVarNames,...
                                'VariableTypes',varTypes,...
                                'Delimiter',delimiter,...
                                'DataLines', dataStartLine,...
                                'ExtraColumnsRule',extraColRule); 

template_log1 = '_photo_results.txt';
template_log2 = '-FID_1_final*.log';
template_rating = '_ratings_results.txt';

nr_runs = 4;
nr_volumes1 = 216; % for run1 and run2
nr_volumes2 = 228; % for run3 and run4

TR = 2.5; % in seconds

% conditions run1 and run2
conditions_name12 = {
    'rest'
    'high calorie'
    'low calorie'
    'neutral'
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

            
%% LOOP OVER SUBJECTS TO READ LOGFILES, CREATE TABLE WITH ONSETS AND DURATIONS, AND SAVE AS EVENTS.TSV TO BIDSSUBJDIRS
%-------------------------------------------------------------------------------------------------------------------
if ~isempty(subjs2write)
    [C,ia,~] = intersect(sourcesubjs,subjs2write);
    
    if ~isequal(C',subjs2write)
        error('\nsubject %s present in subjs2write not present in %s, please check before proceeding',subjs2write{~ismember(subjs2write,C)},derivdir);
    
    else
        
        for sub = ia'
            
            for ses = 1:2
        
            % DEFINE SUBJECT & SESSION LEVEL DIRS
            subjsourcedir = sourcesubjdirs{sub};
            [~,subjectcode,~] = fileparts(subjsourcedir);
            subjBIDSdir = fullfile(BIDSsubjdirs{sub},'func');
            sessubjsourcedir = fullfile(subjsourcedir,['ses-' num2str(ses)]);
            sessubjBIDSdir = fullfile(subjBIDSdir,['ses-' num2str(ses)]);
            logsessubjsourcedir = fullfile(sessubjsourcedir,'logfiles');
            
            go1 = 0; % for run1
            go2 = 0; % for run2
            go3 = 0; % for run3
            go4 = 0; % for run4

                % LOOP OVER RUNS
                for run = 1:size(runnames,2)

                    logfilename = dir(fullfile(logsessubjsourcedir,logfilenames{run}));
                    logfilename = char(logfilename(:).name);
                    logfilepath = fullfile(logsessubjsourcedir,logfilename);

                    if ~isfile(logfilepath)
                        warning('\nlogfile missing for run %d in %s, please check before proceeding',run,logfilepath);
                        continue

                    elseif size(logfilepath,1) > 1
                        error('\nmore than one logfile with run index %s for %s, please check before proceeding',run,sourcesubjs{sub})

                    else
                        
                        if run < 3
                            
                            fid = fopen(logfilepath,'r');
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
                           

                            onset = [];
                            duration = [];
                            trial_type = categorical();
                                for cond = 1:size(conditions_name12,1)
                                    if cond < 5
                                        trial_type_block = categorical();
                                        duration_block = [];
                                        
                                        for block = 1:size(onsets_blocks{cond},1)
                                            trial_type_block(block,1) = conditions_name12{cond};
                                            duration_block(block,1) = durations{cond};
                                        end
                                        
                                        onset = [onset; onsets_blocks{cond}];
                                        duration = [duration; duration_block];
                                        trial_type = [trial_type; trial_type_block];
                                        clear trial_type_block duration_block
                                        
                                    else
                                        trial_type_block = categorical();
                                        
                                        for block = 1:size(onsets_blocks{cond},1)
                                            trial_type_block(block,1) = conditions_name12{cond};
                                        end
                                        
                                        onset = [onset; onsets_blocks{cond}];
                                        duration = [duration; durations{cond}];
                                        trial_type = [trial_type; trial_type_block];
                                        clear trial_type_block
                                        
                                    end
                                end
                                
                             log=table(trial_type,onset,duration);
                             
                                    
                            
                            
                            

                            log.rating = zeros(height(log),1);

                                for l = 1:height(log)

                                    if contains(char(log.Code(l)),'score','IgnoreCase',true)
                                        scorestring = char(log.Code(l));

                                        if strcmp(scorestring(1,end-3:end),'-100')
                                            log.rating(l) = str2double(scorestring(1,end-3:end));

                                        elseif ~contains(scorestring(1,end-2:end),':')
                                            log.rating(l) = str2double(strtrim(scorestring(1,end-2:end)));

                                        else
                                            log.rating(l) = str2double(scorestring(1,end));

                                        end

                                    else
                                        log.rating(l) = NaN;

                                    end

                                end

                            log.trial_type = categorical(log.trial_type);

                            log = log((~isundefined(log.trial_type) | ~isnan(log.rating)),:);

                                for n = 1:height(log)

                                    if ismember(log.trial_type(n),events_interest)
                                        log.rating(n) = log.rating(n+3);

                                    end

                                end

                            log = removevars(log,{'Trial','EventType','Code','Time','TTime','Duration','ReqTime','TimeZero'}); % get rid of junk variables from logfile we don't need

                            log = log(~isundefined(log.trial_type),:);  

                            log.duration = zeros(height(log),1);

                                for m = 1:height(log)

                                    if ~isequal(log.trial_type(m),'fixation')
                                        log.duration(m) = log.onset(m+1) - log.onset(m);

                                    else
                                        log.duration(m) = NaN;

                                    end
                                end

                            log = log(~isnan(log.duration),:);

                            filename = fullfile(subjBIDSdir,[sourcesubjs{sub},'_task-',taskname,runnames{run},'_events.tsv']);
                            writetable(log,filename,'Filetype','text','Delimiter','\t');
                            clear logfile log time_zero filenamee
                            
                        else
                            
                        end

                    end % if loop checking whether logfile exists

                end % for loop runs
            
            end % for loop sessions
            
        end % for loop subjects
    
    end % if loop checking subjs2write present in sourcesubjs        
    
else
    
    if ~isequal(sourcesubjs, BIDSsubjs)
        [D,~,~] = intersect(sourcesubjs,BIDSsubjs);
        error('\nsubject %s present in %s not present in %s, please check before proceeding',BIDSsubjs{~ismember(BIDSsubjs,D)},BIDSdir,derivdir);
        
    else
        
        if pheno_tsv
            pheno_file = table();
            pheno_dir = fullfile(BIDSdir,'phenotype');
                if ~isfolder(pheno_dir)
                    mkdir(pheno_dir);
                end
        end
        
        for sub = 1:size(sourcesubjs,1)
            
            if pheno_tsv
                pheno_file_subj = table();
            end

        % DEFINE SUBJECT LEVEL DIRS
        subjsourcedir = sourcesubjdirs{sub};
        subjBIDSdir = fullfile(BIDSsubjdirs{sub},'func');

            % LOOP OVER RUNS
            for run = 1:size(logfilenames,2)
                
                logfilename = dir(fullfile(subjsourcedir,'logfiles',logfilenames{run}));
                logfilename = char(logfilename(:).name);
                logfilepath = fullfile(subjsourcedir,'logfiles',logfilename);
                
                if ~isfile(logfilepath)
                    warning('\nlogfile missing for run %d in %s, please check before proceeding',run,logfilepath);
                    continue
                
                elseif size(logfilepath,1) > 1
                    error('\nmore than one logfile with run index %s for %s, please check before proceeding',run,sourcesubjs{sub})
                
                else
                    log = readtable(logfilepath,opts);
                    log = log(~isnan(log.Trial),:);
                    time_zero = log.Time(log.Trial == 0 & log.EventType == 'Pulse'); % time for onsets and durations is counted from the first scanner pulse onwards
                        
                        if size(time_zero,1) > 1
                            error('\nambiguity about time zero in %s%s, please check logfile',subjs{sub},logfilenames{run});
                        end
                        
                    log.TimeZero = log.Time - time_zero;
                    log.onset = log.TimeZero ./ 10000; % convert to seconds 
                    log(log.EventType == 'Pulse',:) = [];
                    log.trial_type = cell(height(log),1);
                        
                        for k = 1:height(log)
                            
                            if ismember(log.Code(k),swallow_rinse_labels)
                                log.trial_type{k} = events_nuisance{1};
                            
                            elseif ismember(log.Code(k),rating_labels)
                                log.trial_type{k} = events_nuisance{2};
                                log.onset(k) = log.onset(k)+4;
                            
                            elseif ismember(log.Code(k),sweet_labels)
                                idx = (log.Code(k) == sweet_labels);
                                log.trial_type{k} = events_interest{idx'};
                            
                            elseif ismember(log.Code(k),fixation_labels)
                                log.trial_type{k} = 'fixation';
                            
                            else
                                log.trial_type{k} = '';
                            
                            end
                        end
                    
                    log.rating = zeros(height(log),1);
                    
                        for l = 1:height(log)
                            
                            if contains(char(log.Code(l)),'score','IgnoreCase',true)
                                scorestring = char(log.Code(l));
                                
                                if strcmp(scorestring(1,end-3:end),'-100')
                                    log.rating(l) = str2double(scorestring(1,end-3:end));
                                
                                elseif ~contains(scorestring(1,end-2:end),':')
                                    log.rating(l) = str2double(strtrim(scorestring(1,end-2:end)));
                                
                                else
                                    log.rating(l) = str2double(scorestring(1,end));
                                    
                                end
                            
                            else
                                log.rating(l) = NaN;
                                
                            end
                            
                        end
                        
                    log.trial_type = categorical(log.trial_type);
                    
                    log = log((~isundefined(log.trial_type) | ~isnan(log.rating)),:);
                    
                        for n = 1:height(log)
                            
                            if ismember(log.trial_type(n),events_interest)
                                log.rating(n) = log.rating(n+3);
                            
                            end
                            
                        end
                          
                    log = removevars(log,{'Trial','EventType','Code','Time','TTime','Duration','ReqTime','TimeZero'}); % get rid of junk variables from logfile we don't need
                    
                    log = log(~isundefined(log.trial_type),:);  
                    
                    log.duration = zeros(height(log),1);
                        
                        for m = 1:height(log)
                            
                            if ~isequal(log.trial_type(m),'fixation')
                                log.duration(m) = log.onset(m+1) - log.onset(m);
                            
                            else
                                log.duration(m) = NaN;
                            
                            end
                        end
                        
                    log = log(~isnan(log.duration),:);
                        
                    filename = fullfile(subjBIDSdir,[sourcesubjs{sub},'_task-',taskname,runnames{run},'_events.tsv']);
                    writetable(log,filename,'Filetype','text','Delimiter','\t');
                    
                        if pheno_tsv
                            
                            log2 = log(~isnan(log.rating),:);
                            log2 = removevars(log2,{'onset','duration'});

                            for n = 1:height(log2)
                                log2.participant_id(n,:) = BIDSsubjs{sub};
                                log2.run_id(n) = run;
                                log2.trial_id_run(n) = n; % generates consecutive trial numbers within each run
                                log2.trial_id_concat(n) = height(pheno_file_subj) + n; % generates consecutive trial numbers over all conditions & runs
                            end  
                            
                            for o = 1:size(events_interest,2)
                                idx_run = log2.trial_type == events_interest{o};
                                log2.trial_id_cond_run(idx_run) = 1:sum(idx_run);
                                    if height(pheno_file_subj) > 0
                                        idx_sub = pheno_file_subj.trial_type == events_interest{o};
                                        log2.trial_id_cond_concat(idx_run) = sum(idx_sub)+1:(sum(idx_sub)+sum(idx_run));
                                    else
                                        log2.trial_id_cond_concat = log2.trial_id_cond_run;
                                    end
                                    clear idx_run idx_sub
                            end
                            
                            pheno_file_subj = [pheno_file_subj;log2];

                        end
                    
                    clear logfile log time_zero filename log2

                end % if loop checking whether logfile exists

            end % for loop runs
            
            pheno_file = [pheno_file;pheno_file_subj];
            clear pheno_file_subj;

        end % for loop subjects
        
        pheno_filename = fullfile(pheno_dir,pheno_name);
        writetable(pheno_file,pheno_filename,'Filetype','text','Delimiter','\t');
    
    end % if loop checking sourcesubjs == BIDSsubjs
    
end % if loop checking writing option