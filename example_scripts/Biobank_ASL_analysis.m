% Preliminary analysis of ASL UK Biobank data, Oct 2022

%FigOutDir = '~/scratch/Data/Biobank/ASL_prelim_results/Figs/'

%% Add paths
fsldir = getenv('FSLDIR');
addpath([fsldir '/etc/matlab']);

% Check if we are on FMRIB or BMRC computing clusters
[~,hostname] = system('hostname');
if regexp(hostname,'fmrib') % FMRIB
    localdatdir = '/vols/Data/okell/tokell/Data/Biobank/';
    disp(['FMRIB server detected. Local data dir: ' localdatdir])
elseif regexp(hostname,'bmrc') % BMRC
    localdatdir = '/well/okell/projects/biobank/';
    disp(['BMRC server detected. Local data dir: ' localdatdir])    
else 
    error(['Unknown server location: ' hostname])
end

FigOutDir = [localdatdir 'ASL_analysis/Figs/'];
%stv_dir = '/well/okell/projects/biobank/steve/'; % Local copies of Steve's matlab code from jalapeno, 2/11/22
%addpath([stv_dir 'NETWORKS/FSLNets'],[stv_dir 'matlab'],[stv_dir 'FACS'],[stv_dir 'matlab/FastICA_25']);
% No longer need to add Steve's code as have copied it to my local Matlab
% directory

% Matlab workspace to load variables from:
WSPC = [localdatdir 'workspace13f.mat'];

OUTPUT = [localdatdir 'ASL_analysis/']; % where you want to output stuff

%% Load data
%load([localdatdir 'workspace13f.mat']);
%load([localdatdir 'workspace13f.mat']);
% Skip here and use loadvars below instead to avoid reading everything into
% memory if we don't need it

%% Correlate IDPs with non-imaging measures
% Load the relevant variables
loadvars({'IDPs1_i_deconf','vars_i_deconf','IDP_categories'},WSPC);

% Example code from Steve, June 2022
[unicorrRi_deconf,unicorrPi_deconf,unicorrNi_deconf] = nancorr(IDPs1_i_deconf,vars_i_deconf); % matrix of IDP vs nIDP correlations
grot=-log10(unicorrPi_deconf); grot(IDP_categories~=18,:)=0; % select IDP category for j=1:size(unicorrPi_deconf,2), grot(grot(:,j)<max(grot(:,j)),j)=0; end
 % only keep max IDP association for each nIDP for j=1:size(unicorrPi_deconf,2), for i=1:size(unicorrPi_deconf,1), if varskeepVT(j)>-1 & grot(i,j)>6 % select nIDP categories and logP threshold disp(sprintf('%d %d % .2f %4.1f %4d %s %s',i,j,unicorrRi_deconf(i,j),-log10(unicorrPi_deconf(i,j)),unicorrNi_deconf(i,j),IDP_names{i},varsHeader{j}));
 %disp(sprintf('% .2f %4.1f %4d %s %s %s',unicorrRi_deconf(i,j),-log10(unicorrPi_deconf(i,j)),unicorrNi_deconf(i,j),IDP_names{i},varsVARS{j},varsHeader{j})); end end; end;

 %% Check data in scan 1 or scan 2
 loadvars({'IDPs1_i_deconf','IDPs2_i_deconf','IDP_categories','varsHeader','vars'},WSPC);
 size(IDPs1_i_deconf), size(IDPs2_i_deconf)
 
 ASLScan1Idx = sum(~isnan(IDPs1_i_deconf(:,IDP_categories==18)),2)>0;
 ASLScan2Idx = sum(~isnan(IDPs2_i_deconf(:,IDP_categories==18)),2)>0;
 
 NASLDat1 = sum(ASLScan1Idx)
 NASLDat2 = sum(ASLScan2Idx)
 
 % OK, so ~3500 datasets in scan 1, ~1800 in scan 2
 
 %figure; plot(ASLScan1Idx); hold on; plot(ASLScan2Idx);
 
 % Check for overlap
 NASLScan1AndScan2 = sum(ASLScan1Idx & ASLScan2Idx)
 
 % Look at Covid status of scan 1 and scan 2 data
 grotCv=unique([nets_cellfind(varsHeader,'covid',1)]);
 %varsHeader(grotCv)'  % Sanity check 
 
 grotCvCaseControl = nets_cellfind(varsHeader,'Case-control status for COVID19 imaging repeat (3.0)');
 %grotCvVirusIdent = nets_cellfind(varsHeader,'Diagnoses - ICD10 (U071 - U07.1 COVID-19 virus identified)');
 %grotCvVirusNotIdent = nets_cellfind(varsHeader,'Diagnoses - ICD10 (U072 - U07.2 COVID-19 virus not identified)');
 varsHeader(grotCvCaseControl)'  % Sanity check 
 %varsHeader(grotCvVirusIdent)'  % Sanity check 
 %varsHeader(grotCvVirusNotIdent)'  % Sanity check 
 
% OK, these don't seem to tally up as expected - CaseControl gives an
% integer - turns out this is coded - the final digit tells you whether
% it's a control (0) or Covid case (1).
CV = mod(vars(:,grotCvCaseControl),10);

disp(['Number of ASL scan 1 data sets with Covid control = ' ns(sum(CV==0 & ASLScan1Idx))])
disp(['Number of ASL scan 1 data sets with Covid cases = '   ns(sum(CV==1 & ASLScan1Idx))])

disp(['Number of ASL scan 2 data sets with Covid control = ' ns(sum(CV==0 & ASLScan2Idx))])
disp(['Number of ASL scan 2 data sets with Covid cases = '   ns(sum(CV==1 & ASLScan2Idx))])

% Virus identified and not identified is binary, but very few not
% identifieds, so presumably relates to virus subtypes being identifiable
% or not, or similar.

% Could consider merging scan1 and scan2 (non-Covid) cases, but wouldn't
% change numbers dramatically: 3414 vs. 4287, so perhaps leave for now.

 %% ASL results
 % Dimensions of grot are IDPs x nIDPs
 % Number of statistical tests for ASL is no ASL IDPs x no nIDPs
 NumASLIDPs = sum(IDP_categories == 18);
 NumnIDPs = size(vars_i_deconf,2);
 NumTest = NumASLIDPs * NumnIDPs;
 BonfThr = -log10(0.05/NumTest);
 
 % FDR correction
 pVals = unicorrPi_deconf(IDP_categories==18,:); % Grab ASL p values
 pVals = pVals(~isnan(pVals));
 [pID pN] = FDR(pVals(:),0.05);
 FDRThr = -log10(pN); % Take the more conservative assumption to use as a threshold
 
 %% Look at only the max IDP associations for each nIDP
 grot2 = grot;
 for j=1:size(unicorrPi_deconf,2), grot2(grot2(:,j)<max(grot2(:,j)),j)=0; end % only keep max IDP association for each nIDP 
 for j=1:size(unicorrPi_deconf,2), 
     for i=1:size(unicorrPi_deconf,1), 
         if varskeepVT(j)>-1 & grot2(i,j)>FDRThr % select nIDP categories and logP threshold 
             %disp(sprintf('%d %d % .2f %4.1f %4d %s %s',i,j,unicorrRi_deconf(i,j),-log10(unicorrPi_deconf(i,j)),unicorrNi_deconf(i,j),IDP_names{i},varsHeader{j}));
             disp(sprintf('% .2f %4.1f %4d %s %s %s',       unicorrRi_deconf(i,j),-log10(unicorrPi_deconf(i,j)),unicorrNi_deconf(i,j),IDP_names{i},varsVARS{j},varsHeader{j})); 
         end
     end
 end
 
 %% Extract nIDP category names from Matlab script
 %txt = fileread('/vols/Data/ukbiobank/FMRIB/ANALYSIS/DO_13/subproc_2_nIDP_groups.m');
 txt = fileread([localdatdir 'DO_13/subproc_2_nIDP_groups.m']);
 % Discard up to phenotype categories
 txt = regexprep(txt,'.*phenotype categories','');
 % Discard beyond categories
 txt = regexprep(txt,'%%%% lists of used and not used categories.*','');
 
 disp(txt)
 nIDPcatnames = {};
 for ii = 1:100
     ln = regexp(txt,[' ' ns(ii) ' [^\n]*\n'],'match');
     nIDPcatnames{ii} = regexprep(regexprep(ln,[' ' ns(ii) ' *'],''),'\n*','');
     %nIDPcatnames{ii} = regexprep(regexprep(txt,['.*% ' ns(ii) ' *'],''),'\n.*','');
 end
 
 %% Plot by nIDP category
 xstart = 1; leg = {}; jj = 1;
 LargeFigWindow(1,1); cols = distinguishable_colors(100,[1 1 1]);
 for ii = 1:max(varskeepVT)
     % Grab the relevant results
     y = grot(IDP_categories==18,varskeepVT==ii);
     
     % Remove those nIDPs where there are no results > 0 or non-NaN
     y = y(:,sum(y,1)>0);
     
     n = size(y,2);
     
     if n > 0
         x = xstart:(xstart+n-1);
         %y = grot(IDP_categories==18,varskeepVT==ii);
         
         x = repmat(x,size(y,1),1);
         scatter(x(:),y(:),5,cols(jj,:),'filled'); hold on;
         xstart = max(x(:))+1;
         leg{jj} = nIDPcatnames{ii}{1};
         jj = jj+1;
     end
 end


xlabel 'nIDPs'
ylabel '-log10(p)';
xlim([0 max(x(:))])
plot(xlim,[FDRThr FDRThr],'k--');
plot(xlim,[BonfThr BonfThr],'k:');
leg{end+1} = 'FDR';
leg{end+1} = 'Bonferroni';
legend(leg);

toexportfig(gcf,'ASL_Manhatten_plot',FigOutDir,1);

%% Deconfound with haematocrit
% Find haematocrit related nIDPs
%grotHct=unique([nets_cellfind(varsHeader,'haematocrit',1) nets_cellfind(varsHeader,'haemoglobin',1)]);

%% Look for other haemoglobin type measures in case there are more of these
grotHaem=unique([nets_cellfind(varsHeader,'Haematocrit',1)]);

varsHeader(grotHaem)' % List outputs

for ii = grotHaem
    disp(varsHeader(ii))
    disp([sum(~isnan(vars(:,ii)))/length(vars(:,ii)) ...
    sum(~isnan(vars(:,ii)) & ASLScan1Idx)/sum(ASLScan1Idx) ...
    sum(~isnan(vars(:,ii)) & ASLScan2Idx)/sum(ASLScan2Idx) ])
end

% OK, most subjects have 'Haematocrit percentage (0.0)'

% ??Perhaps modify to use the most recent Haematocrit measurement? Unclear
% why so few subjects had follow-up measurements...

%% Look at haematocrit values and correlation with GM CBF
loadvars({'varsHeader','vars','IDP_names','IDPs1'},WSPC);

HctIdx=unique([nets_cellfind(varsHeader,'Haematocrit percentage (0.0)',1)]);

varsHeader(HctIdx)'  % Sanity check

figure; histogram(vars(:,HctIdx));

% Look at Haemotrocrit vs. mean GM CBF
GMCBFIdx = nets_cellfind(IDP_names,'ASL_region_analysis_gm_-_70%_GM');
IDP_names{GMCBFIdx}


x = vars(:,HctIdx);
y = IDPs1(:,GMCBFIdx);
Idx = ~isnan(y);
figure; plot(x(Idx),y(Idx),'.');
xlabel 'Haematocrit'
ylabel 'Mean GM CBF'

% Negative correlation, as expected.

%% Try deconfounding with haematocrit
loadvars({'conf_names','conf1','subject_IDs_unique','subject_IDs','site'},WSPC);

% Extract haematocrit
confHct_raw = vars(:,HctIdx);

% Extract site for scan one, using code similar to DO_1_load_IDPs.m
[~,grotI1,grotJ1]=intersect(subject_IDs_unique+2e7,subject_IDs);

site1=zeros(size(subject_IDs_unique,1),1)/0; 
site1(grotI1,:)=site(grotJ1,:); 

% Create cell array of subject indices by site
subjectIndicesBySite = {};
for SiteNo = 1:max(site1(:))
    subjectIndicesBySite{SiteNo} = find(site1 == SiteNo);
end

% Demedian, split by site and remove outliers, as per Fidel's deconfounding
% paper
[confHct, confHct_names] = toduplicateDemedianNormBySite({'Hct_perc'},confHct_raw,subject_IDs_unique, subjectIndicesBySite);

% Check mean is zero and std is one (for non-zero entries)
disp(['Mean confHct: ' ns(mean(confHct))])
for ii = 1:3; tmp = confHct(:,ii); tmpstd(ii) = std(tmp(tmp~=0)); end
disp(['Std confHct: ' ns(tmpstd)])

% Look at the histogram before and after demedianing etc.
figure; subplot(1,3,1); 
histogram(confHct_raw); title 'Hct percentage (raw)';
% NaNs don't show up here - if we set to the median we can see them
subplot(1,3,2);
disp([ns(sum(isnan(confHct_raw))) ' NaNs in the original Hct data'])
tmp = confHct_raw; tmp(isnan(tmp)) = median(tmp(~isnan(tmp)));
histogram(tmp); title 'Hct percentage (NaN set to median)';
subplot(1,3,3); 
histogram(confHct(confHct~=0)); title 'Hct percentage (normalised)';
% Can also plot by site:
%for SiteNo = 1:max(site1(:))
%    tmp = confHct(:,SiteNo); 
%    histogram(tmp(subjectIndicesBySite{SiteNo})); hold on;
%end


%% Construct the confound matrix
% Ideally should run this on all the data, but just stick to scan 1 for
% simplicity here
conf1_haem = [conf1 confHct];
conf1_haem_names = [conf_names; mat2cell(confHct_names,ones(size(confHct_names,1),1),size(confHct_names,2))];

%% Deconfound only the ASL IDPs
loadvars({'IDPs1_i_deconf','IDP_categories','IDPs1_i'},WSPC);
IDPs1_i_deconf_haem=IDPs1_i_deconf;
IDPs1_i_deconf_haem(:,IDP_categories==18) = nets_unconfound_par(IDPs1_i(:,IDP_categories==18),conf1_haem);

%% Repeat analysis above
loadvars('vars_i_deconf',WSPC);
% Set Hct related variables to NaN to avoid spurious correlations after
% deconfounding
HctAllIdx = unique([nets_cellfind(varsHeader,'Haematocrit',1)]);
varsHeader(HctAllIdx)
vars_i_deconf_haem = vars_i_deconf; 
vars_i_deconf_haem(:,HctAllIdx) = NaN;
[unicorrRi_deconf_haem,unicorrPi_deconf_haem,unicorrNi_deconf_haem] = nancorr(IDPs1_i_deconf_haem,vars_i_deconf_haem); % matrix of IDP vs nIDP correlations
grot_haem=-log10(unicorrPi_deconf_haem); grot_haem(IDP_categories~=18,:)=0; % select IDP category for j=1:size(unicorrPi_deconf,2), grot(grot(:,j)<max(grot(:,j)),j)=0; end

%% Recalculate Bonf and FDR thresholds
NumASLIDPs = sum(IDP_categories == 18);
 NumnIDPs = size(vars_i_deconf,2);
 NumTest = NumASLIDPs * NumnIDPs;
 BonfThr_haem = -log10(0.05/NumTest);
 
 % FDR correction
 pVals_haem = unicorrPi_deconf_haem(IDP_categories==18,:); % Grab ASL p values
 pVals_haem= pVals_haem(~isnan(pVals_haem));
 [pID_haem pN_haem] = FDR(pVals_haem(:),0.05);
 FDRThr_haem = -log10(pN_haem); % Take the more conservative assumption to use as a threshold
 

%% Plot by nIDP category
loadvars('varskeepVT',WSPC);
 xstart = 1; leg = {}; jj = 1;
 LargeFigWindow(1,1); cols = distinguishable_colors(100,[1 1 1]);
 for ii = 1:max(varskeepVT)
     % Grab the relevant results
     y = grot_haem(IDP_categories==18,varskeepVT==ii);
     
     % Remove those nIDPs where there are no results > 0 or non-NaN
     y = y(:,sum(y,1)>0);
     
     n = size(y,2);
     
     if n > 0
         x = xstart:(xstart+n-1);
         %y = grot(IDP_categories==18,varskeepVT==ii);
         
         x = repmat(x,size(y,1),1);
         scatter(x(:),y(:),5,cols(jj,:),'filled'); hold on;
         xstart = max(x(:))+1;
         leg{jj} = nIDPcatnames{ii}{1};
         jj = jj+1;
     end
 end


xlabel 'nIDPs'
ylabel '-log10(p)';
xlim([0 max(x(:))])
plot(xlim,[FDRThr_haem FDRThr_haem],'k--');
plot(xlim,[BonfThr_haem BonfThr_haem],'k:');
leg{end+1} = 'FDR';
leg{end+1} = 'Bonferroni';
legend(leg);

toexportfig(gcf,'ASL_Manhatten_plot_haematocrit_deconfounded',FigOutDir,1);







%% ----- Deconfound with haematocrit and cortical thickness ----- %%
% Look for cortical thickness IDPs
grotCortthk=unique([nets_cellfind(IDP_names,'GlobalMeanThickness',-1)]);

IDP_names(grotCortthk) % List outputs

%% Look at left and right histograms
figure; histogram(IDPs1(:,grotCortthk(1))); hold on;
histogram(IDPs1(:,grotCortthk(2))); 
legend(IDP_names(grotCortthk))

%% Take the mean
mean_cortthk = mean(IDPs1(:,grotCortthk),2);

%% Check correlation with age
loadvars('age1',WSPC);
figure; scatter(age1,mean_cortthk,'.');

% Small negative correlation

%% Demedian, remove outliers etc.
[confcortthk, confcortthk_names] = toduplicateDemedianNormBySite({'MeanCorticalThickness'},mean_cortthk,subject_IDs_unique, subjectIndicesBySite);

% Check mean is zero and std is one (for non-zero entries)
disp(['Mean confcortthk: ' ns(mean(confcortthk))])
for ii = 1:3; tmp = confcortthk(:,ii); tmpstd(ii) = std(tmp(tmp~=0)); end
disp(['Std confcortthk: ' ns(tmpstd)])

% Look at the histogram before and after demedianing etc.
figure; subplot(1,3,1); 
histogram(mean_cortthk); title 'cortthk (raw)';
% NaNs don't show up here - if we set to the median we can see them
subplot(1,3,2);
disp([ns(sum(isnan(mean_cortthk))) ' NaNs in the original cortthk data'])
tmp = mean_cortthk; tmp(isnan(tmp)) = median(tmp(~isnan(tmp)));
histogram(tmp); title 'cortthk (NaN set to median)';
subplot(1,3,3); 
histogram(confcortthk(confcortthk~=0)); title 'cortthk (normalised)';
% Can also plot by site:
%for SiteNo = 1:max(site1(:))
%    tmp = confcortthk(:,SiteNo); 
%    histogram(tmp(subjectIndicesBySite{SiteNo})); hold on;
%end


%% Construct the confound matrix
% Ideally should run this on all the data, but just stick to scan 1 for
% simplicity here
conf1_haem_cortthk = [conf1_haem confcortthk];
conf1_haem_cortthk_names = [conf1_haem_names; mat2cell(confcortthk_names,ones(size(confcortthk_names,1),1),size(confcortthk_names,2))];

%% Deconfound only the ASL IDPs
IDPs1_i_deconf_haem_cortthk=IDPs1_i_deconf;
IDPs1_i_deconf_haem_cortthk(:,IDP_categories==18) = nets_unconfound_par(IDPs1_i(:,IDP_categories==18),conf1_haem_cortthk);

%% Repeat analysis above
% Remove Hct from the vars list first
HctAllIdx = unique([nets_cellfind(varsHeader,'Haematocrit',1)]);
varsHeader(HctAllIdx)
vars_i_deconf_haem_cortthk = vars_i_deconf; 
vars_i_deconf_haem_cortthk(:,HctAllIdx) = NaN;
[unicorrRi_deconf_haem_cortthk,unicorrPi_deconf_haem_cortthk,unicorrNi_deconf_haem_cortthk] = nancorr(IDPs1_i_deconf_haem_cortthk,vars_i_deconf_haem_cortthk); % matrix of IDP vs nIDP correlations
grot_haem_cortthk=-log10(unicorrPi_deconf_haem_cortthk); grot_haem_cortthk(IDP_categories~=18,:)=0; % select IDP category for j=1:size(unicorrPi_deconf,2), grot(grot(:,j)<max(grot(:,j)),j)=0; end

%% Exclude some categories to match Miller et al. Nat Neuro 2016
ExcList = {}; %{'genetics','hearing test','eye test','health and medical history','mental health self-report'};
ExcnIDPs = 0;
for ii = 1:length(ExcList)
    ExcIdx = nets_cellfind(nIDPcatnames,ExcList{ii});
    if length(ExcIdx) ~= 1
        error(['Not only one variable for ' ExcList{ii}])
    end
    grot_haem_cortthk(:,varskeepVT==ExcIdx) = 0;
    ExcnIDPs = ExcnIDPs + sum(varskeepVT==ExcIdx);
end

%% Recalculate Bonf and FDR thresholds
NumASLIDPs = sum(IDP_categories == 18);
 NumnIDPs = size(vars_i_deconf,2) - ExcnIDPs;
 NumTest = NumASLIDPs * NumnIDPs;
 BonfThr_haem_cortthk = -log10(0.05/NumTest);
 
 % FDR correction
 pVals_haem_cortthk = unicorrPi_deconf_haem_cortthk(IDP_categories==18,:); % Grab ASL p values
 pVals_haem_cortthk= pVals_haem_cortthk(~isnan(pVals_haem_cortthk));
 [pID_haem_cortthk pN_haem_cortthk] = FDR(pVals_haem_cortthk(:),0.05);
 FDRThr_haem_cortthk = -log10(pN_haem_cortthk); % Take the more conservative assumption to use as a threshold
 
 %% Clear unnecessary variables
 clear vars_i_deconf_haem_cortthk
 
 %% Save
 save([OUTPUT 'ASL_nIDP_associations_haematocrit_and_cort_thk_deconfounded.mat'],'*_haem_cortthk*')

 %% Look at only the max IDP associations for each nIDP
 loadvars('varsVARS',WSPC);
 fileID = fopen(['ASL_nIDP_max_IDP_associations_haematocrit_and_cort_thk_deconfounded.txt'],'w');
 grot2 = grot_haem_cortthk;
 for j=1:size(unicorrPi_deconf_haem_cortthk,2), grot2(grot2(:,j)<max(grot2(:,j)),j)=0; end % only keep max IDP association for each nIDP 
 for j=1:size(unicorrPi_deconf_haem_cortthk,2), 
     for i=1:size(unicorrPi_deconf_haem_cortthk,1), 
         if varskeepVT(j)>-1 & grot2(i,j)>FDRThr_haem_cortthk % select nIDP categories and logP threshold 
             %disp(sprintf('%d %d % .2f %4.1f %4d %s %s',i,j,unicorrRi_deconf(i,j),-log10(unicorrPi_deconf(i,j)),unicorrNi_deconf(i,j),IDP_names{i},varsHeader{j}));
             disp(sprintf('% .2f %4.1f %4d %s %s %s',       unicorrRi_deconf_haem_cortthk(i,j),-log10(unicorrPi_deconf_haem_cortthk(i,j)),unicorrNi_deconf_haem_cortthk(i,j),IDP_names{i},varsVARS{j},varsHeader{j})); 
             fprintf(fileID,'% .2f %4.1f %4d %s %s %s \n',       unicorrRi_deconf_haem_cortthk(i,j),-log10(unicorrPi_deconf_haem_cortthk(i,j)),unicorrNi_deconf_haem_cortthk(i,j),IDP_names{i},varsVARS{j},varsHeader{j});
         end
     end
 end

fclose(fileID);
 
%% Plot by nIDP category
 xstart = 1; leg = {}; jj = 1;
 LargeFigWindow(1,1); cols = distinguishable_colors(100,[1 1 1]);
 for ii = 1:max(varskeepVT)
     % Grab the relevant results
     y = grot_haem_cortthk(IDP_categories==18,varskeepVT==ii);
     
     % Remove those nIDPs where there are no results > 0 or non-NaN
     y = y(:,sum(y,1)>0);
     
     n = size(y,2);
     
     if n > 0
         x = xstart:(xstart+n-1);
         %y = grot(IDP_categories==18,varskeepVT==ii);
         
         x = repmat(x,size(y,1),1);
         %scatter(x(:),y(:),5,cols(jj,:),'filled'); hold on;
         scatter(x(:),y(:),25,cols(jj,:)); hold on;
         xstart = max(x(:))+1;
         leg{jj} = nIDPcatnames{ii}{1};
         jj = jj+1;
     end
 end


xlabel 'nIDPs'
ylabel '-log10(p)';
xlim([0 max(x(:))])
plot(xlim,[FDRThr_haem_cortthk FDRThr_haem_cortthk],'k--');
plot(xlim,[BonfThr_haem_cortthk BonfThr_haem_cortthk],'k:');
leg{end+1} = 'FDR';
leg{end+1} = 'Bonferroni';
legend(leg,'location','eastoutside');

toexportfig(gcf,'ASL_Manhatten_plot_haematocrit_and_cort_thk_deconfounded',FigOutDir,0.5);







%% ------ Look at correlations with other IDPs ------- %%
[unicorrRi_deconfIDP,unicorrPi_deconfIDP,unicorrNi_deconfIDP] = nancorr(IDPs1_i_deconf,IDPs1_i_deconf); % matrix of IDP vs IDP correlations
grotIDP=-log10(unicorrPi_deconfIDP); grotIDP(IDP_categories~=18,:)=0; % select IDP category for j=1:size(unicorrPi_deconf,2), grot(grot(:,j)<max(grot(:,j)),j)=0; end
% Remove self-correlations
grotIDP(IDP_categories==18,IDP_categories==18) = 0;

%% Recalculate Bonferroni and FDR thresholds
% Dimensions of grot are IDPs x IDPs
 % Number of statistical tests for ASL is no ASL IDPs x no non-ASL IDPs
 NumASLIDPs = sum(IDP_categories == 18);
 NumnonASLIDPs = sum(IDP_categories ~=18);
 NumTestIDP = NumASLIDPs * NumnonASLIDPs;
 BonfThrIDP = -log10(0.05/NumTestIDP);
 
 % FDR correction
 pValsIDP = unicorrPi_deconfIDP(IDP_categories==18,:); % Grab ASL p values
 pValsIDP = pValsIDP(~isnan(pValsIDP));
 [pIDIDP pNIDP] = FDR(pValsIDP(:),0.05);
 FDRThrIDP = -log10(pNIDP); % Take the more conservative assumption to use as a threshold
 
  %% Look at only the max non-ASL IDP associations for each ASL IDP
 grot2IDP = grotIDP;
 for j=1:size(unicorrPi_deconfIDP,2), grot2IDP(grot2IDP(:,j)<max(grot2IDP(:,j)),j)=0; end % only keep max ASL IDP association for each non-ASL IDP 
 for j=1:size(unicorrPi_deconfIDP,2), 
     for i=1:size(unicorrPi_deconfIDP,1), 
         if varskeepVT(j)>-1 & grot2IDP(i,j)>FDRThrIDP % select nIDP categories and logP threshold 
             %disp(sprintf('%d %d % .2f %4.1f %4d %s %s',i,j,unicorrRi_deconf(i,j),-log10(unicorrPi_deconf(i,j)),unicorrNi_deconf(i,j),IDP_names{i},varsHeader{j}));
             disp(sprintf('% .2f %4.1f %4d %s %s ',       unicorrRi_deconfIDP(i,j),-log10(unicorrPi_deconfIDP(i,j)),unicorrNi_deconfIDP(i,j),IDP_names{i},IDP_names{j})); 
         end
     end
 end
 
%% Plot
xstart = 1; leg = {}; jj = 1; Ncats = max(IDP_categories(:));
 LargeFigWindow(1,1); cols = distinguishable_colors(Ncats,[1 1 1]);
 for ii = 1:Ncats
     % Grab the relevant results
     y = grotIDP(IDP_categories==18,IDP_categories==ii);
     
     % Remove those IDPs where there are no results > 0 or non-NaN
     y = y(:,sum(y,1)>0);
     
     n = size(y,2);
     
     if n > 0
         x = xstart:(xstart+n-1);
         %y = grot(IDP_categories==18,varskeepVT==ii);
         
         x = repmat(x,size(y,1),1);
         %scatter(x(:),y(:),5,cols(jj,:),'filled'); hold on;
         scatter(x(:),y(:),25,cols(jj,:)); hold on;
         xstart = max(x(:))+1;
         leg{jj} = IDP_category_names{ii};
         jj = jj+1;
     end
 end


xlabel 'IDPs'
ylabel '-log10(p)';
xlim([0 max(x(:))])
plot(xlim,[FDRThrIDP FDRThrIDP],'k--');
plot(xlim,[BonfThrIDP BonfThrIDP],'k:');
leg{end+1} = 'FDR';
leg{end+1} = 'Bonferroni';
legend(leg,'location','eastoutside');

toexportfig(gcf,'ASL_Manhatten_plot_IDPs',FigOutDir,0.5);



%% Look at correlations with other IDPs after deconfounding with haematocrit and cortical thickness
??Remove cortical thickness before running this!!
[unicorrRi_deconfIDP_haem_cortthk,unicorrPi_deconfIDP_haem_cortthk,unicorrNi_deconfIDP_haem_cortthk] = nancorr(IDPs1_i_deconf_haem_cortthk,IDPs1_i_deconf_haem_cortthk); % matrix of IDP vs IDP correlations
grotIDP_haem_cortthk=-log10(unicorrPi_deconfIDP_haem_cortthk); grotIDP_haem_cortthk(IDP_categories~=18,:)=0; % select IDP category for j=1:size(unicorrPi_deconf,2), grot(grot(:,j)<max(grot(:,j)),j)=0; end
% Remove self-correlations
grotIDP_haem_cortthk(IDP_categories==18,IDP_categories==18) = 0;

%% Recalculate Bonferroni and FDR thresholds
% Dimensions of grot are IDPs x IDPs
 % Number of statistical tests for ASL is no ASL IDPs x no non-ASL IDPs
 NumASLIDPs = sum(IDP_categories == 18);
 NumnonASLIDPs = sum(IDP_categories ~=18);
 NumTestIDP = NumASLIDPs * NumnonASLIDPs;
 BonfThrIDP_haem_cortthk = -log10(0.05/NumTestIDP);
 
 % FDR correction
 pValsIDP_haem_cortthk = unicorrPi_deconfIDP_haem_cortthk(IDP_categories==18,:); % Grab ASL p values
 pValsIDP_haem_cortthk = pValsIDP_haem_cortthk(~isnan(pValsIDP_haem_cortthk));
 [pIDIDP_haem_cortthk pNIDP_haem_cortthk] = FDR(pValsIDP_haem_cortthk(:),0.05);
 FDRThrIDP_haem_cortthk = -log10(pNIDP_haem_cortthk); % Take the more conservative assumption to use as a threshold
 
  %% Look at only the max non-ASL IDP associations for each ASL IDP
 fileID = fopen(['ASL_IDP_max_IDP_associations_haematocrit_and_cort_thk_deconfounded.txt'],'w');
 grot2IDP_haem_cortthk = grotIDP_haem_cortthk;
 for j=1:size(unicorrPi_deconfIDP_haem_cortthk,2), grot2IDP_haem_cortthk(grot2IDP_haem_cortthk(:,j)<max(grot2IDP_haem_cortthk(:,j)),j)=0; end % only keep max ASL IDP association for each non-ASL IDP 
 for j=1:size(unicorrPi_deconfIDP_haem_cortthk,2), 
     for i=1:size(unicorrPi_deconfIDP_haem_cortthk,1), 
         if varskeepVT(j)>-1 & grot2IDP_haem_cortthk(i,j)>FDRThrIDP_haem_cortthk % select nIDP categories and logP threshold 
             %disp(sprintf('%d %d % .2f %4.1f %4d %s %s',i,j,unicorrRi_deconf(i,j),-log10(unicorrPi_deconf(i,j)),unicorrNi_deconf(i,j),IDP_names{i},varsHeader{j}));
             disp(sprintf('% .2f %4.1f %4d %s %s ',       unicorrRi_deconfIDP_haem_cortthk(i,j),-log10(unicorrPi_deconfIDP_haem_cortthk(i,j)),unicorrNi_deconfIDP_haem_cortthk(i,j),IDP_names{i},IDP_names{j})); 
             fprintf(fileID,'% .2f %4.1f %4d %s %s \n',unicorrRi_deconfIDP_haem_cortthk(i,j),-log10(unicorrPi_deconfIDP_haem_cortthk(i,j)),unicorrNi_deconfIDP_haem_cortthk(i,j),IDP_names{i},IDP_names{j});
         end
     end
 end

fclose(fileID);

%% Plot
xstart = 1; leg = {}; jj = 1; Ncats = max(IDP_categories(:));
 LargeFigWindow(1,1); cols = distinguishable_colors(Ncats,[1 1 1]);
 for ii = 1:Ncats
     % Grab the relevant results
     y = grotIDP_haem_cortthk(IDP_categories==18,IDP_categories==ii);
     
     % Remove those IDPs where there are no results > 0 or non-NaN
     y = y(:,sum(y,1)>0);
     
     n = size(y,2);
     
     if n > 0
         x = xstart:(xstart+n-1);
         %y = grot(IDP_categories==18,varskeepVT==ii);
         
         x = repmat(x,size(y,1),1);
         %scatter(x(:),y(:),5,cols(jj,:),'filled'); hold on;
         scatter(x(:),y(:),25,cols(jj,:)); hold on;
         xstart = max(x(:))+1;
         leg{jj} = IDP_category_names{ii};
         jj = jj+1;
     end
 end


xlabel 'IDPs'
ylabel '-log10(p)';
xlim([0 max(x(:))])
plot(xlim,[FDRThrIDP_haem_cortthk FDRThrIDP_haem_cortthk],'k--');
plot(xlim,[BonfThrIDP_haem_cortthk BonfThrIDP_haem_cortthk],'k:');
leg{end+1} = 'FDR';
leg{end+1} = 'Bonferroni';
legend(leg,'location','eastoutside');

toexportfig(gcf,'ASL_Manhatten_plot_IDPs_haematocrit_and_cort_thk_deconfounded',FigOutDir,0.5);






%% ----- Look at correlation with age ------- %%
% Need to first deconfound with haematocrit and cortical thickness without
% age
grot=nets_cellfind(conf_names,'Age',-1); conf_names(grot)  % find all age-related confounds

IDPs1_i_deconf_haem_cortthk_noage=IDPs1_i_deconfnoage;
IDPs1_i_deconf_haem_cortthk_noage(:,IDP_categories==18) = nets_unconfound_par(IDPs1_i(:,IDP_categories==18),conf1_haem_cortthk(:, setdiff([1:size(conf1_haem_cortthk,2)],grot)));

%% Try only with haematocrit
grot=nets_cellfind(conf_names,'Age',-1); conf_names(grot)  % find all age-related confounds

IDPs1_i_deconf_haem_noage=IDPs1_i_deconfnoage;
IDPs1_i_deconf_haem_noage(:,IDP_categories==18) = nets_unconfound_par(IDPs1_i(:,IDP_categories==18),conf1_haem(:, setdiff([1:size(conf1_haem,2)],grot)));


%% Plot for each ASL IDP
ASL_IDP_Idx = find(IDP_categories==18);
cols = distinguishable_colors(length(ASL_IDP_Idx),[1 1 1]);
figure;
for ii = 1:1; %length(ASL_IDP_Idx)
    IDPNo = ASL_IDP_Idx(ii);
    %Idx = ~isnan(age1(:)) & ~isnan(IDPs1_i_deconf_haem_cortthk_noage(:,IDPNo));
    %Idx = ~isnan(age1(:)) & ~isnan(IDPs1_i_deconf_haem_noage(:,IDPNo));
    %dx = ~isnan(age1(:)) & ~isnan(IDPs1_i_deconfnoage(:,IDPNo));
    Idx = ~isnan(age1(:)) & ~isnan(IDPs1(:,IDPNo));
    x = age1(Idx);
    %y = IDPs1_i_deconf_haem_cortthk_noage(Idx,IDPNo); % Minimal correlation
    %y = IDPs1_i_deconf_haem_noage(Idx,IDPNo); % Minimal correlation
    %y = IDPs1_i_deconfnoage(Idx,IDPNo); % Minimal correlation
    %y = IDPs1_i(Idx,IDPNo); % This still shows the correlation
    y = IDPs1(Idx,IDPNo);
    %scatter(x,y,5,cols(ii,:),'filled'); hold on;
    [xp, ymean, ysd] = to_sliding_mean_sd(x,y,5,100);
    shadedErrorBar(xp,ymean,ysd,'lineProps',{'markerfacecolor',cols(ii,:)}); hold on;
end
legend(IDP_names(ASL_IDP_Idx));

%% Just mean GM CBF
IDPNo = nets_cellfind(IDP_names,'ASL_region_analysis_gm_-_70%_GM');
figure;

Idx = ~isnan(age1(:)) & ~isnan(IDPs1(:,IDPNo));
x = age1(Idx);
y = IDPs1(Idx,IDPNo);
[xp, ymean, ysd] = to_sliding_mean_sd(x,y,5,100);
shadedErrorBar(xp,ymean,ysd,'lineProps',{'markerfacecolor',cols(ii,:),'linewidth',2}); hold on;
xlabel 'Age/years'; ylabel 'Mean GM CBF (ml/100g/min)'

toexportfig(gcf,'CBF_vs_age_plot',FigOutDir,0.2)

%% Just mean GM arrival
IDPNo = nets_cellfind(IDP_names,'ASL_region_analysis_arrival_gm_-_70%_GM');
figure;

Idx = ~isnan(age1(:)) & ~isnan(IDPs1(:,IDPNo));
x = age1(Idx);
y = IDPs1(Idx,IDPNo);
[xp, ymean, ysd] = to_sliding_mean_sd(x,y,5,100);
shadedErrorBar(xp,ymean,ysd,'lineProps',{'color',cols(2,:),'linewidth',2}); hold on;
xlabel 'Age/years'; ylabel 'Mean GM ATT (s)'

toexportfig(gcf,'ATT_vs_age_plot',FigOutDir,0.2)


%% Consider why trend is removed after deconfounding
% It seems unconfounding with cortical thickness removes almost all age
% related ASL signal... In fact, even without cortical thickness or
% haematocrit there is almost no age-related change in mean GM CBF...
% Look for correlations between age and other confounding factors
grot=nets_cellfind(conf_names,'Age',-1); conf_names(grot)  % find all age-related confounds
Idx = setdiff([1:size(conf1,2)],grot);
conf1_noage = conf1(:, Idx);
conf1_noage_names = conf_names(Idx);

[Rage,Page,Nage] = nancorr(conf1_noage,age1);

% TODO:
% Calculate moving average and SD
% Average R and L regions
% Use GM Only
% Plot ATT separately






%% ------ Generate spatial maps for interesting associations -----%%
% Clear memory since we will need it below!
%toclear;

% Reload variables we need for now
%load([localdatdir 'workspace13f.mat'],'subject_IDs_unique','IDPs1_i_deconf','IDP_categories');
loadvars({'subject_IDs_unique','IDPs1_i_deconf','IDP_categories'},WSPC);
ASLScan1Idx = sum(~isnan(IDPs1_i_deconf(:,IDP_categories==18)),2)>0;
clear IDPs1_i_deconf

% Based on example code from Chaoyue
MODALITY='ASL';
MASK=[getenv('FSLDIR') '/data/standard/MNI152_T1_2mm_brain_mask.nii.gz'];
UKBDATA = '/well/win-biobank/projects/imaging/data/data3/subjectsAll/';

%% Check data is present for subjects that have some non-NaN ASL IDPs
subj_IDs_ASL = subject_IDs_unique(ASLScan1Idx);
subjs_with_ASL = 0;
for ii = 1:length(subj_IDs_ASL)
    if ~exist([UKBDATA '2' ns(subj_IDs_ASL(ii)) '/ASL/BASIL/OXASL_ra_dir/std_space/perfusion_calib.nii.gz'],'file')
   %if ~exist([UKBDATA '2' ns(subj_IDs_ASL(ii)) '/ASL/BASIL/OXASL_ra_dir/std_space/arrival.nii.gz'],'file')
        disp(['Data does not exist for subject ' ns(subj_IDs_ASL(ii))])
    else
        subjs_with_ASL = subjs_with_ASL + 1;
    end
end

disp(['Subjects with ASL data = ' ns(subjs_with_ASL) ' out of ' ns(length(subj_IDs_ASL)) ' expected'])

% All seems to be present

%% Read in the mask
%%addpath /well/win/projects/ukbiobank/fbp/bb_FSL/etc/matlab
%clear tempdir; setenv('TMPDIR','/dev/shm'); tempdir;

%cd(UKBDATA)
%grot=sprintf('ls -1 */%s.nii.gz | awk -F / ''{print $1}'' | sort -n > %s/%s.txt',MODALITY,OUTPUT,MODALITY)
%system(grot)
%cd(OUTPUT)
mask=niftiread(MASK);  maskORIG=mask;  mask=mask>0;  mask=reshape(mask,prod(size(mask)),1);
%IDs=load(sprintf('%s.txt',MODALITY));

%% Check one subject
testdat = single(niftiread([datdir '2' ns(subj_IDs_ASL(1)) '/ASL/BASIL/OXASL_ra_dir/std_space/perfusion_calib.nii.gz']));
DispIm([single(maskORIG(:,:,round(end/2))) testdat(:,:,round(end/2))/max(testdat(:))])
clear testdat
% Looks good

%% Read in all the data and save out
N=length(subj_IDs_ASL);
BIG=zeros(N,sum(mask),'single')/0;
%FTYPE='perfusion_calib';
FTYPE='arrival';

for i=1:N
  disp([i N]);
  FName = [datdir '2' ns(subj_IDs_ASL(i)) '/ASL/BASIL/OXASL_ra_dir/std_space/' FTYPE '.nii.gz']
  if exist(FName)>0
    grot=niftiread(FName);  BIG(i,:)=grot(mask>0)';
  end
end

clear grot*; 

save([OUTPUT FTYPE '_AllSubjs4D.mat'],'-v7.3');

%% Reload if needed
%FTYPE = 'perfusion_calib';
FTYPE = 'arrival';
load([OUTPUT FTYPE '_AllSubjs4D.mat'])

%% Regenerate the confounding factors with haematocrit and cortical thickness
% Reload variables we need for now
%load([localdatdir 'workspace13f.mat'],'varsHeader','vars','conf1','conf_names','IDP_names','IDPs1');
loadvars({'varsHeader','vars','conf1','conf_names','IDP_names','IDPs1'},WSPC);

HctIdx=unique([nets_cellfind(varsHeader,'Haematocrit percentage (0.0)',1)]);
conf1_haem = [conf1 vars(:,HctIdx)];
conf1_haem_names = conf_names;
conf1_haem_names{end+1} = 'Haematocrit percentage (0.0)';
grotCortthk=unique([nets_cellfind(IDP_names,'GlobalMeanThickness',-1)]);
mean_cortthk = mean(IDPs1(:,grotCortthk),2);
conf1_haem_cortthk = [conf1_haem mean_cortthk];
conf1_haem_cortthk_names = conf1_haem_names;
conf1_haem_cortthk_names{end+1} = 'Mean cortical thickness';
clear varsHeader vars conf1 conf_names IDP_names IDPs1 conf1_haem conf1_haem_names

%% Save out
save([OUTPUT 'conf1_haem_cortthk.mat'],'conf1_haem_cortthk','conf1_haem_cortthk_names','ASLScan1Idx','-v7.3');

%% Reload if needed
load([OUTPUT 'conf1_haem_cortthk.mat']);

%% Extract confounding factors for the relevant subjects
conf=conf1_haem_cortthk(ASLScan1Idx,:);
clear conf1_haem_cortthk

%% Quantile normalisation
disp('Running quantile normalisation!');
MyPhenos = nets_inormal( BIG );
save([OUTPUT FTYPE '_AllSubjs4D_i.mat'],'MyPhenos','-v7.3');
clear BIG

%% Reload if needed
load([OUTPUT FTYPE '_AllSubjs4D_i.mat']);

%% Deconfound
disp('Running deconfounding...')
MyPhenos_Deconf = to_nets_unconfound2(MyPhenos,conf);
save([OUTPUT FTYPE '_AllSubjs4D_i_deconf.mat'],'MyPhenos_Deconf','-v7.3');
clear MyPhenos
% This takes a very long time... Try splitting into blocks

%% Test one block
ASL_spatial_deconfound_block(1,2000,'perfusion_calib')

%% Generate scripts
cd(OUTPUT)
!NBlocks=1000; for (( ii=1; ii<=$NBlocks; ii++ )); do echo "addpath ~/Documents/Matlab/; run ~/Documents/Matlab/startup.m; ASL_spatial_deconfound_block($ii,$NBlocks,'perfusion_calib')" > blocks/block${ii}of${NBlocks}.m; done

%% Try array approach for arrival
!rm blocks/arrival_array.txt
!NBlocks=1000; for (( ii=1; ii<=$NBlocks; ii++ )); do echo "matlab -nojvm -nodisplay -nosplash -singleCompThread -r \"addpath ~/Documents/Matlab/; run ~/Documents/Matlab/startup.m; ASL_spatial_deconfound_block($ii,$NBlocks,'arrival')\" " >> blocks/arrival_array.txt; done

%% Some of these failed - re-run manually
% FJOBS=$(for (( ii=1; ii<=1000; ii++ )); do if [ ! -e blocks/perfusion_calib_AllSubjs4D_i_deconf_block${ii}of1000.mat ]; then echo $ii; fi; done)
% for ii in $FJOBS; do BatchMatlabJob short.qc blocks/block${ii}of1000.m; done
% FJOBS=$(for (( ii=1; ii<=1000; ii++ )); do if [ ! -e blocks/arrival_AllSubjs4D_i_deconf_block${ii}of1000.mat ]; then echo $ii; fi; done)
% for ii in $FJOBS; do fsl_sub -q short.qc $(sed -n ${ii}p blocks/arrival_array.txt ); done

%% Check how long it takes to calculate the correlation for one block
%load([localdatdir 'workspace13f.mat'],'vars_i_deconf','varsHeader');
loadvars({'vars_i_deconf','varsHeader'},WSPC);

nIDPIdx=unique([nets_cellfind(varsHeader,'Cardiac output during PWA (2.0)',1)]);
varsHeader(nIDPIdx)

niIDPDeconf=vars_i_deconf(ASLScan1Idx,nIDPIdx); % take only one nIDP

tmp = load('blocks/perfusion_calib_AllSubjs4D_i_deconf_block500of1000.mat');

tic
[unicorrRi_deconf_block,unicorrPi_deconf_block,unicorrNi_deconf_block] = nancorr(tmp.MyPhenos_Deconf_block,niIDPDeconf);
toc
% Very fast - let's reconstitute the blocks of data first

%% Reconstitute the blocks of data
NBlocks = 1000;
FTYPE = 'arrival';
MyPhenos_Deconf = zeros(size(MyPhenos))/0; % start with NaN values
startcol = 1;
for ii = 1:NBlocks
    disp(ii)
    tmp = load(['blocks/' FTYPE '_AllSubjs4D_i_deconf_block' ns(ii) 'of' ns(NBlocks) '.mat']);
    Ncols = size(tmp.MyPhenos_Deconf_block,2);
    MyPhenos_Deconf(:,startcol:(startcol+Ncols-1)) = tmp.MyPhenos_Deconf_block;
    startcol = startcol + Ncols;
end

% Save
save([OUTPUT FTYPE '_AllSubjs4D_i_deconf.mat'],'MyPhenos_Deconf','-v7.3');

%% Reload if needed
FTYPE = 'arrival';
load([OUTPUT FTYPE '_AllSubjs4D_i_deconf.mat']);

%% Find the nIDPs with strongest associations with ASL IDPs
if ~exist('grot_haem_cortthk','var')
    load([localdatdir 'ASL_nIDP_associations_haematocrit_and_cort_thk_deconfounded.mat'],'grot_haem_cortthk');
end

 grot2 = grot_haem_cortthk;
 nIDPs_of_interest_Idx = []; nIDPs_of_interest_names = {};
 for j=1:size(grot_haem_cortthk,2), grot2(grot2(:,j)<max(grot2(:,j)),j)=0; end % only keep max IDP association for each nIDP 
 for j=1:size(grot_haem_cortthk,2), 
     for i=1:size(grot_haem_cortthk,1), 
         if varskeepVT(j)>-1 & grot2(i,j)>FDRThr_haem_cortthk % select nIDP categories and logP threshold 
             %disp(sprintf('%d %d % .2f %4.1f %4d %s %s',i,j,unicorrRi_deconf(i,j),-log10(unicorrPi_deconf(i,j)),unicorrNi_deconf(i,j),IDP_names{i},varsHeader{j}));
             nIDPs_of_interest_Idx = [nIDPs_of_interest_Idx j];
             nIDPs_of_interest_names{end+1} = varsHeader{j};              
         end
     end
 end
clear grot2
varsHeader(nIDPs_of_interest_Idx)

%% Write the nIDPs to a file
fileID = fopen([FigOutDir 'ASL_nIDPs_of_interest.txt'],'w');
for ii = 1:length(nIDPs_of_interest_Idx)
    fprintf(fileID,'%d %s \n', ii-1, varsHeader{nIDPs_of_interest_Idx(ii)});
end

fclose(fileID);
%% Run the correlation for all these IDPs of interest
% if ~exist('vars_i_deconf','var')
%     load([localdatdir 'workspace13f.mat'],'vars_i_deconf','varsHeader');
% end
loadvars({'vars_i_deconf','varsHeader'},WSPC);

%nIDPIdx=unique([nets_cellfind(varsHeader,'Cardiac output during PWA (2.0)',1)]);
%nIDP_name = varsHeader(nIDPIdx)

%niIDPDeconf=vars_i_deconf(ASLScan1Idx,nIDPIdx); % take only one nIDP
niIDPDeconf=vars_i_deconf(ASLScan1Idx,nIDPs_of_interest_Idx);

tic
[unicorrRi_deconf_haem_cortthk_vox,unicorrPi_deconf_haem_cortthk_vox,unicorrNi_deconf_haem_cortthk_vox] = nancorr(MyPhenos_Deconf,niIDPDeconf);
toc

%% Convert back to image space
unicorrRi_deconf_haem_cortthk_vox_map = zeros([size(maskORIG) length(nIDPs_of_interest_Idx)]);
unicorrPi_deconf_haem_cortthk_vox_map = zeros([size(maskORIG) length(nIDPs_of_interest_Idx)]);
unicorrRi_deconf_haem_cortthk_vox_map(repmat(maskORIG > 0,[1 1 1 length(nIDPs_of_interest_Idx)])) = unicorrRi_deconf_haem_cortthk_vox(:);
unicorrPi_deconf_haem_cortthk_vox_map(repmat(maskORIG > 0,[1 1 1 length(nIDPs_of_interest_Idx)])) = unicorrPi_deconf_haem_cortthk_vox(:);

%% Save out
%OUTNAME = [FTYPE '_' regexprep(nIDP_name{1},'[ ().]','_') '_R'];
OUTNAME = [FTYPE '_All_nIDPs_of_interest_R'];
save_avw(unicorrRi_deconf_haem_cortthk_vox_map,OUTNAME,'f',[1 1 1 1]);

tosystem(['fslcpgeom ' MASK ' ' OUTNAME ' -d'])

%% View
tosystem(['fsleyes -std ' OUTNAME ' -dr 0.03 0.1 -un -cm hot -nc blue-lightblue'])





%% -------- Repeat for IDPs ------------ %%
% Find the IDPs with strongest associations with ASL IDPs
if ~exist('grotIDP_haem_cortthk','var')
    error('regenerate grotIDP_haem_cortthk');
end

 grot2IDP_haem_cortthk = grotIDP_haem_cortthk;
 IDPs_of_interest_Idx = []; IDPs_of_interest_names = {};
 for j=1:size(grotIDP_haem_cortthk,2), grot2IDP_haem_cortthk(grot2IDP_haem_cortthk(:,j)<max(grot2IDP_haem_cortthk(:,j)),j)=0; end % only keep max IDP association for each ASL IDP 
 for j=1:size(grotIDP_haem_cortthk,2), 
     for i=1:size(grotIDP_haem_cortthk,1), 
         if grot2IDP_haem_cortthk(i,j)>FDRThrIDP_haem_cortthk % select nIDP categories and logP threshold 
             %disp(sprintf('%d %d % .2f %4.1f %4d %s %s',i,j,unicorrRi_deconf(i,j),-log10(unicorrPi_deconf(i,j)),unicorrNi_deconf(i,j),IDP_names{i},varsHeader{j}));
             IDPs_of_interest_Idx = [IDPs_of_interest_Idx j];
             IDPs_of_interest_names{end+1} = IDP_names{j};              
         end
     end
 end
%clear grot2IDP_haem_cortthk
IDPs_of_interest_names'

%% Only keep the strongest association in each category to narrow it down a bit
IDPs_of_interest_cat_Idx = [];
IDPs_of_interest_cat_names = {};
jj = 1;
for ii = unique(IDP_categories')
    maxgrot = max(max(grot2IDP_haem_cortthk(:,IDP_categories==ii)));
    if maxgrot > FDRThrIDP_haem_cortthk
        IDPs_of_interest_cat_Idx(jj) = find((max(grot2IDP_haem_cortthk)==maxgrot) & (IDP_categories'==ii));
        IDPs_of_interest_cat_names{end+1} = IDP_names{IDPs_of_interest_cat_Idx(jj)};
        jj = jj + 1;
    end
end

IDPs_of_interest_cat_names'

%% Write the IDPs to a file
fileID = fopen([FigOutDir 'ASL_IDPs_of_interest.txt'],'w');
for ii = 1:length(IDPs_of_interest_cat_Idx)
    fprintf(fileID,'%d %s \n', ii-1, IDP_names{IDPs_of_interest_cat_Idx(ii)});
end

fclose(fileID);
%% Run the correlation for all these IDPs of interest
if ~exist('IDPs1_i_deconf_haem_cortthk','var')
    error('IDPs1_i_deconf_haem_cortthk not present...')
end

%FTYPE = 'perfusion_calib';
FTYPE = 'arrival';
load([OUTPUT FTYPE '_AllSubjs4D_i_deconf.mat'],'MyPhenos_Deconf');

iIDPDeconf=IDPs1_i_deconf_haem_cortthk(ASLScan1Idx,IDPs_of_interest_cat_Idx);

tic
[unicorrRi_deconf_haem_cortthk_vox,unicorrPi_deconf_haem_cortthk_vox,unicorrNi_deconf_haem_cortthk_vox] = nancorr(MyPhenos_Deconf,iIDPDeconf);
toc

%% Convert back to image space
unicorrRi_deconf_haem_cortthk_vox_map = zeros([size(maskORIG) length(IDPs_of_interest_cat_Idx)]);
unicorrPi_deconf_haem_cortthk_vox_map = zeros([size(maskORIG) length(IDPs_of_interest_cat_Idx)]);
unicorrRi_deconf_haem_cortthk_vox_map(repmat(maskORIG > 0,[1 1 1 length(IDPs_of_interest_cat_Idx)])) = unicorrRi_deconf_haem_cortthk_vox(:);
unicorrPi_deconf_haem_cortthk_vox_map(repmat(maskORIG > 0,[1 1 1 length(IDPs_of_interest_cat_Idx)])) = unicorrPi_deconf_haem_cortthk_vox(:);

%% Save out
%OUTNAME = [FTYPE '_' regexprep(nIDP_name{1},'[ ().]','_') '_R'];
OUTNAME = [FTYPE '_All_IDPs_of_interest_R'];
save_avw(unicorrRi_deconf_haem_cortthk_vox_map,OUTNAME,'f',[1 1 1 1]);

tosystem(['fslcpgeom ' MASK ' ' OUTNAME ' -d'])

%% View
tosystem(['fsleyes -std ' OUTNAME ' -dr 0.03 0.1 -un -cm hot -nc blue-lightblue'])






%% ----- Look at correlations with age -----%%
% Reload data
FTYPE = 'perfusion_calib';
load([OUTPUT FTYPE '_AllSubjs4D_i.mat']);

%% Set up confounds
grot=nets_cellfind(conf1_haem_names,'Age',-1); conf1_haem_names(grot)  % find all age-related confounds
% Make sure haematocrit percentage is not included
Idx = nets_cellfind(conf1_haem_names(grot),'Haematocrit percentage')
grot = grot(setdiff([1:length(grot)],Idx)); conf1_haem_names(grot)

conf1_haem_noage = conf1_haem(:, setdiff([1:size(conf1_haem,2)],grot));
conf1_haem_noage_names = conf1_haem_names(setdiff([1:size(conf1_haem,2)],grot));

%% Save
save([OUTPUT 'conf1_haem_noage.mat'],'conf1_haem_noage','conf1_haem_noage_names','ASLScan1Idx','-v7.3');

%% Reload if needed
load([OUTPUT 'conf1_haem_noage.mat']);

%% Run deconfounding
% Extract the relevant confounding data
conftmp=conf1_haem_noage(ASLScan1Idx,:);

disp('Running deconfounding...')
MyPhenos_Deconf = to_nets_unconfound2(MyPhenos,conftmp);

save([OUTPUT FTYPE '_AllSubjs4D_i_deconf_noage.mat'],'MyPhenos_Deconf','-v7.3');
clear MyPhenos
% Still very slow - split into blocks

%% Test one block
ASL_spatial_deconfound_noage_block(1,1000,'perfusion_calib')

%% Generate scripts
cd(OUTPUT)
!rm blocks/perfusion_calib_noage_array.txt
!NBlocks=1000; for (( ii=1; ii<=$NBlocks; ii++ )); do echo "matlab -nojvm -nodisplay -nosplash -singleCompThread -r \"addpath ~/Documents/Matlab/; run ~/Documents/Matlab/startup.m; ASL_spatial_deconfound_noage_block($ii,$NBlocks,'perfusion_calib')\" " >> blocks/perfusion_calib_noage_array.txt; done
% Submit with: fsl_sub -q short.qc -t ./blocks/perfusion_calib_noage_array.txt

%% Repeat for arrival
!rm blocks/arrival_noage_array.txt
!NBlocks=1000; for (( ii=1; ii<=$NBlocks; ii++ )); do echo "matlab -nojvm -nodisplay -nosplash -singleCompThread -r \"addpath ~/Documents/Matlab/; run ~/Documents/Matlab/startup.m; ASL_spatial_deconfound_noage_block($ii,$NBlocks,'arrival')\" " >> blocks/arrival_noage_array.txt; done
% Submit with: fsl_sub -q short.qc -t ./blocks/arrival_noage_array.txt

%% Reconstitute the blocks of data
NBlocks = 1000;
FTYPE = 'arrival';
MyPhenos_Deconf_noage = zeros(size(MyPhenos))/0; % start with NaN values
startcol = 1;
for ii = 1:NBlocks
    disp(ii)
    tmp = load(['blocks/' FTYPE '_AllSubjs4D_i_deconf_noage_block' ns(ii) 'of' ns(NBlocks) '.mat']);
    Ncols = size(tmp.MyPhenos_Deconf_noage_block,2);
    MyPhenos_Deconf_noage(:,startcol:(startcol+Ncols-1)) = tmp.MyPhenos_Deconf_noage_block;
    startcol = startcol + Ncols;
end

% Save
save([OUTPUT FTYPE '_AllSubjs4D_i_deconf_noage.mat'],'MyPhenos_Deconf_noage','-v7.3');

%% Reload if needed
%FTYPE = 'perfusion_calib';
FTYPE = 'arrival';
%load([OUTPUT FTYPE '_AllSubjs4D_i_deconf_noage.mat']);
load([OUTPUT FTYPE '_AllSubjs4D_i.mat']);


%% Run the correlation with age
% NB. use non-deconfounded data since unconfounding even with non-age
% variables still removes main correlations with age because there are
% correlations between age and other confounding factors it seems

% if ~exist('age1','var')
%     load([localdatdir 'workspace13f.mat'],'age1');
% end
loadvars('age1',WSPC);

tic
%[unicorrRi_deconf_haem_noage_vox,unicorrPi_deconf_haem_noage_vox,unicorrNi_deconf_haem_noage_vox] = nancorr(MyPhenos_Deconf_noage,age1(ASLScan1Idx));
[unicorrRi_age_vox,unicorrPi_age_vox,unicorrNi_age_vox] = nancorr(MyPhenos,age1(ASLScan1Idx));
toc

%% Convert back to image space
unicorrRi_age_vox_map = zeros([size(maskORIG) 1]);
unicorrPi_age_vox_map = zeros([size(maskORIG) 1]);
unicorrRi_age_vox_map(repmat(maskORIG > 0,[1 1 1 1])) = unicorrRi_age_vox(:);
unicorrPi_age_vox_map(repmat(maskORIG > 0,[1 1 1 1])) = unicorrPi_age_vox(:);

%% Save out
%OUTNAME = [FTYPE '_' regexprep(nIDP_name{1},'[ ().]','_') '_R'];
OUTNAME = [FTYPE '_age_R'];
save_avw(unicorrRi_age_vox_map,OUTNAME,'f',[1 1 1 1]);

tosystem(['fslcpgeom ' MASK ' ' OUTNAME ' -d'])

%% View
tosystem(['fsleyes -std ' OUTNAME ' -dr 0.03 0.1 -un -cm hot -nc blue-lightblue'])


%% Also generate mean and std maps across subjects
% Reload data
FTYPE = 'arrival';
load([FTYPE '_AllSubjs4D.mat']);

%% Take the mean and SD ignoring NaNs
MeanIm = nanmean(BIG,1);
SDIm   = nanstd(BIG,[],1);

%% Save out
MeanIm_map = zeros([size(maskORIG) 1]);
SDIm_map   = zeros([size(maskORIG) 1]);
MeanIm_map(maskORIG > 0) = MeanIm(:);
SDIm_map(  maskORIG > 0) = SDIm(:);

%% Save out
%OUTNAME = [FTYPE '_' regexprep(nIDP_name{1},'[ ().]','_') '_R'];
OUTNAME = [FTYPE '_mean'];
save_avw(MeanIm_map,OUTNAME,'f',[1 1 1 1]);
tosystem(['fslcpgeom ' MASK ' ' OUTNAME ' -d'])

OUTNAME = [FTYPE '_SD'];
save_avw(SDIm_map,OUTNAME,'f',[1 1 1 1]);
tosystem(['fslcpgeom ' MASK ' ' OUTNAME ' -d'])



%% -----Covid vs. non-Covid, out of interest----- %%
% Loop through ASL IDPs
AllIDPs = 1:length(IDP_names);
ASLIDPs = AllIDPs(IDP_categories==18);
CovTStatsP = zeros(length(ASLIDPs),1);
CovTStatsP_nodeconf = zeros(length(ASLIDPs),1);
figure;
for ii = 1:length(ASLIDPs)
    IDPNo = ASLIDPs(ii);
    ControlIDPs = IDPs2_i_deconf(ASLScan2Idx & (CV==0),IDPNo);
    CaseIDPs = IDPs2_i_deconf(ASLScan2Idx & (CV==1),IDPNo);
    [~,CovTStatsP(ii)] = ttest2(ControlIDPs,CaseIDPs,'Vartype','unequal');
    ControlIDPs_nodeconf = IDPs2(ASLScan2Idx & (CV==0),IDPNo);
    ControlIDPs_nodeconf = ControlIDPs_nodeconf(~isnan(ControlIDPs_nodeconf));
    CaseIDPs_nodeconf = IDPs2(ASLScan2Idx & (CV==1),IDPNo);
    CaseIDPs_nodeconf = CaseIDPs_nodeconf(~isnan(CaseIDPs_nodeconf));
    [~,CovTStatsP_nodeconf(ii)] = ttest2(ControlIDPs_nodeconf,CaseIDPs_nodeconf,'Vartype','unequal');
%     errorbar(ii-0.2,mean(ControlIDPs(~isnan(ControlIDPs))),std(ControlIDPs(~isnan(ControlIDPs))),'b.'); hold on;
%     errorbar(ii+0.2,mean(CaseIDPs(~isnan(CaseIDPs))),std(CaseIDPs(~isnan(CaseIDPs))),'r.'); hold on;
    if regexpi(IDP_names{IDPNo},'arrival')
        subplot(1,2,2); title 'ATT';
    else % CBF
        subplot(1,2,1); title 'CBF';
    end
    errorbar(ii-0.2,mean(ControlIDPs_nodeconf),std(ControlIDPs_nodeconf),'b.'); hold on;
    errorbar(ii+0.2,mean(CaseIDPs_nodeconf),std(CaseIDPs_nodeconf),'r.'); hold on;
end

% Definite tendency for CBF to be lower and ATT to be longer in Covid
% subjects, but small compared to variance
BonfThrCov = -log10(0.05/length(ASLIDPs));
 
 % FDR correction
[pIDCov, pNCov] = FDR(CovTStatsP(:),0.05);
FDRThrCov = -log10(pNCov); % Take the more conservative assumption to use as a threshold
 
figure; plot(-log10(CovTStatsP)); hold on;
plot([1 length(ASLIDPs)],[BonfThrCov BonfThrCov],'k--');
title 'Deconfounded'

% Repeat for not deconfounded
 % FDR correction
[pIDCov_nodeconf, pNCov_nodeconf] = FDR(CovTStatsP_nodeconf(:),0.05);
FDRThrCov_nodeconf = -log10(pNCov_nodeconf); % Take the more conservative assumption to use as a threshold
 
figure; plot(-log10(CovTStatsP_nodeconf)); hold on;
plot([1 length(ASLIDPs)],[BonfThrCov BonfThrCov],'k--');
title 'Not deconfounded'

% Nothing very significant here, esp. after multiple comparisons correction



%% Voxelwise comparison of Covid vs. non-Covid
% Used FEAT to generate example .con and .mat files in:
% /well/okell/projects/biobank/ASL_analysis/Covid/Example_FEAT_setup

%% Create 4D data
loadvars('subject_IDs_unique',WSPC);

% Find subject IDs
CovCntlIDs = subject_IDs_unique(ASLScan2Idx & (CV==0));
CovCaseIDs = subject_IDs_unique(ASLScan2Idx & (CV==1));
AllCovIDs = [CovCntlIDs(:); CovCaseIDs(:)];
AllCovIdx = ASLScan2Idx & ~isnan(CV);

% Based on example code from Chaoyue
MODALITY='ASL';
MASK=[getenv('FSLDIR') '/data/standard/MNI152_T1_2mm_brain_mask.nii.gz'];
UKBDATA = '/well/win-biobank/projects/imaging/data/data3/subjectsAll/';

%% Check data is present for subjects that have some non-NaN ASL IDPs
subjs_with_ASL = 0;
for ii = 1:length(AllCovIDs)
    if ~exist([UKBDATA '3' ns(AllCovIDs(ii)) '/ASL/BASIL/OXASL_ra_dir/std_space/perfusion_calib.nii.gz'],'file')
   %if ~exist([UKBDATA '2' ns(subj_IDs_ASL(ii)) '/ASL/BASIL/OXASL_ra_dir/std_space/arrival.nii.gz'],'file')
        disp(['Data does not exist for subject ' ns(AllCovIDs(ii))])
    else
        subjs_with_ASL = subjs_with_ASL + 1;
    end
end

disp(['Subjects with ASL data = ' ns(subjs_with_ASL) ' out of ' ns(length(AllCovIDs)) ' expected'])

% All seems to be present

%% Read in the mask
mask=niftiread(MASK);  maskORIG=mask;  mask=mask>0;  mask=reshape(mask,prod(size(mask)),1);

%% Check one subject
testdat = single(niftiread([UKBDATA '3' ns(AllCovIDs(1)) '/ASL/BASIL/OXASL_ra_dir/std_space/perfusion_calib.nii.gz']));
DispIm([single(maskORIG(:,:,round(end/2))) testdat(:,:,round(end/2))/max(testdat(:))])
clear testdat
% Looks good

%% Read in all the data and save out
N=length(AllCovIDs);
BIG=zeros(N,sum(mask),'single')/0;
%FTYPE='perfusion_calib';
FTYPE='arrival';

for i=1:N
  disp([i N]);
  FName = [UKBDATA '3' ns(AllCovIDs(i)) '/ASL/BASIL/OXASL_ra_dir/std_space/' FTYPE '.nii.gz']
  if exist(FName)>0
    grot=niftiread(FName);  BIG(i,:)=grot(mask>0)';
  end
end

clear grot*; 

save([OUTPUT 'Covid/' FTYPE '_AllCovSubjs4D.mat'],'BIG','-v7.3');

%% Reload if needed
%FTYPE = 'perfusion_calib';
FTYPE = 'arrival';
load([OUTPUT 'Covid/' FTYPE '_AllCovSubjs4D.mat'])

%% Save out as a 4D file for randomise
% Convert back to 4D data
BIG4D = zeros([size(maskORIG) length(AllCovIDs)]);
% Note BIG has dimensions subjects by voxels, so we need to permute to get
% voxels first, then subjects
BIG4D(repmat(maskORIG > 0,[1 1 1 length(AllCovIDs)])) = permute(BIG,[2 1]);

% Save
OUTNAME = [OUTPUT 'Covid/' FTYPE '_AllCovSubjs4D'];
save_avw(BIG4D,OUTNAME,'f',[1 1 1 1]);
tosystem(['fslcpgeom ' MASK ' ' OUTNAME ' -d'])

%% Set up randomise design files
OUTDIR = [OUTPUT 'Covid/' FTYPE '_Cov_randomise'];
if ~exist(OUTDIR); mkdir(OUTDIR); end

% Copy the contrasts file from the example
tosystem(['cp ' OUTPUT 'Covid/Example_FEAT_setup/cov_vs_con.con ' OUTDIR '/design.con'])

% Write the design.mat file
DESF = [OUTDIR '/design.mat'];
fid = fopen(DESF,'w');
fprintf(fid,'%s\n','/NumWaves       2');
fprintf(fid,'%s\n',['/NumPoints      ' ns(length(AllCovIDs))]);
fprintf(fid,'%s\n','/PPheights              1.000000e+00    1.000000e+00');
fprintf(fid,'%s\n','');
fprintf(fid,'%s\n','/Matrix');
for ii=1:length(AllCovIDs)
    % Find the subject index
    Idx = find(subject_IDs_unique == AllCovIDs(ii));
    
    % Write format: Covid (0=no, 1=yes)   Control (0=no, 1=yes)
    fprintf(fid,'%d\t%d\n',CV(Idx),~CV(Idx));
end
fclose(fid);

%% Run randomise with a small number of permutations to check it runs
tosystem(['randomise -i ' OUTPUT 'Covid/' FTYPE '_AllCovSubjs4D -o ' OUTDIR '/cov_vs_con_quickT -d ' OUTDIR '/design.mat -t ' OUTDIR '/design.con -m ' getenv('FSLDIR') '/data/standard/MNI152_T1_2mm_brain_mask.nii.gz -n 50 -T'])

%% Full run
tosystem(['fsl_sub -q long.q randomise -i ' OUTPUT 'Covid/' FTYPE '_AllCovSubjs4D -o ' OUTDIR '/cov_vs_con_' FTYPE 'T -d ' OUTDIR '/design.mat -t ' OUTDIR '/design.con -m ' getenv('FSLDIR') '/data/standard/MNI152_T1_2mm_brain_mask.nii.gz -n 1000 -T'])

% Some spatial patterns, esp. reduced CBF near left frontal/olfactory
% regions, but not significant.

%% Regenerate the confounding factors with haematocrit and cortical thickness for scan 2
% Reload variables we need for now
%load([localdatdir 'workspace13f.mat'],'varsHeader','vars','conf2','conf_names','IDP_names','IDPs2');
loadvars({'varsHeader','vars','conf2','conf_names','IDP_names','IDPs2'},WSPC);

HctIdx=unique([nets_cellfind(varsHeader,'Haematocrit percentage (0.0)',1)]);
conf2_haem = [conf2 vars(:,HctIdx)];
conf2_haem_names = conf_names;
conf2_haem_names{end+1} = 'Haematocrit percentage (0.0)';
grotCortthk=unique([nets_cellfind(IDP_names,'GlobalMeanThickness',-1)]);
mean_cortthk = mean(IDPs2(:,grotCortthk),2); %??USE IDPs2_i??
conf2_haem_cortthk = [conf2_haem mean_cortthk];
conf2_haem_cortthk_names = conf2_haem_names;
conf2_haem_cortthk_names{end+1} = 'Mean cortical thickness';
clear varsHeader vars conf2 conf_names IDP_names IDPs1 conf2_haem conf2_haem_names

%% Save out
save([OUTPUT 'conf2_haem_cortthk.mat'],'conf2_haem_cortthk','conf2_haem_cortthk_names','ASLScan2Idx','AllCovIdx','-v7.3');

%% Reload if needed
load([OUTPUT 'conf2_haem_cortthk.mat']);

%% Extract confounding factors for the Covid subjects
conf=conf2_haem_cortthk(AllCovIdx,:);
clear conf2_haem_cortthk

%% Quantile normalisation
disp('Running quantile normalisation!');
MyPhenos = nets_inormal( BIG );
save([OUTPUT 'Covid/' FTYPE '_AllCovSubjs4D_i.mat'],'MyPhenos','-v7.3');
clear BIG

%% Reload if needed
load([OUTPUT 'Covid/' FTYPE '_AllCovSubjs4D_i.mat']);

%% Deconfound
disp('Running deconfounding...')
MyPhenos_Deconf = to_nets_unconfound2(MyPhenos,conf);
save([OUTPUT 'Covid/' FTYPE '_AllCovSubjs4D_i_deconf.mat'],'MyPhenos_Deconf','-v7.3');
clear MyPhenos
% This takes a very long time... Try splitting into blocks

%% Test one block
ASL_spatial_deconfound_block(1,2000,'perfusion_calib',true)

%% Generate scripts
%cd(OUTPUT)
%!NBlocks=1000; for (( ii=1; ii<=$NBlocks; ii++ )); do echo "addpath ~/Documents/Matlab/; run ~/Documents/Matlab/startup.m; ASL_spatial_deconfound_block($ii,$NBlocks,'perfusion_calib')" > blocks/block${ii}of${NBlocks}.m; done

%% Try array approach
cd([OUTPUT '/Covid'])
tosystem(['rm blocks/' FTYPE '_array.txt']);
%tosystem(['NBlocks=1000; for (( ii=1; ii<=$NBlocks; ii++ )); do echo "matlab -nojvm -nodisplay -nosplash -singleCompThread -r \"addpath ~/Documents/Matlab/; run ~/Documents/Matlab/startup.m; ASL_spatial_deconfound_block($ii,$NBlocks,''' FTYPE ''')\" " >> blocks/' FTYPE '_array.txt; done']);
NBlocks = 1000;
fid = fopen(['blocks/' FTYPE '_array.txt'],'w');
for ii = 1:NBlocks
    fprintf(fid,"%s\n",['matlab -nojvm -nodisplay -nosplash -singleCompThread -r "addpath ~/Documents/Matlab/; run ~/Documents/Matlab/startup.m; ASL_spatial_deconfound_block(' ns(ii) ',' ns(NBlocks) ',''' FTYPE ''',true)" ']);
end
fclose(fid);

%% Run
% fsl_sub -q veryshort.q -t blocks/perfusion_calib_array.txt 
% fsl_sub -q veryshort.q -t blocks/arrival_array.txt 

%% Reconstitute the blocks of data
%FTYPE = 'perfusion_calib';
FTYPE = 'arrival';
MyPhenos_Deconf = zeros(size(MyPhenos))/0; % start with NaN values
startcol = 1;
for ii = 1:NBlocks
    disp(ii)
    tmp = load(['blocks/' FTYPE '_AllCovSubjs4D_i_deconf_block' ns(ii) 'of' ns(NBlocks) '.mat']);
    Ncols = size(tmp.MyPhenos_Deconf_block,2);
    MyPhenos_Deconf(:,startcol:(startcol+Ncols-1)) = tmp.MyPhenos_Deconf_block;
    startcol = startcol + Ncols;
end

% Save
save([OUTPUT 'Covid/' FTYPE '_AllCovSubjs4D_i_deconf.mat'],'MyPhenos_Deconf','-v7.3');

%% Reload if needed
FTYPE = 'arrival';
load([OUTPUT 'Covid/' FTYPE '_AllCovSubjs4D_i_deconf.mat']);

%% Save as Nifti
% Convert back to 4D data
MyPhenos_Deconf4D = zeros([size(maskORIG) length(AllCovIDs)]);
% Note MyPhenos_Deconf4D has dimensions subjects by voxels, so we need to permute to get
% voxels first, then subjects
MyPhenos_Deconf4D(repmat(maskORIG > 0,[1 1 1 length(AllCovIDs)])) = permute(MyPhenos_Deconf,[2 1]);

% Save
OUTNAME = [OUTPUT 'Covid/' FTYPE '_AllCovSubjs4D_i_deconf'];
save_avw(MyPhenos_Deconf4D,OUTNAME,'f',[1 1 1 1]);
tosystem(['fslcpgeom ' MASK ' ' OUTNAME ' -d'])

%% Run randomise on this deconfounded data
OUTDIR = [OUTPUT 'Covid/' FTYPE '_Cov_randomise'];
tosystem(['fsl_sub -q long.q randomise -i ' OUTPUT 'Covid/' FTYPE '_AllCovSubjs4D_i_deconf -o ' OUTDIR '/cov_vs_con_' FTYPE '_i_deconfT -d ' OUTDIR '/design.mat -t ' OUTDIR '/design.con -m ' getenv('FSLDIR') '/data/standard/MNI152_T1_2mm_brain_mask.nii.gz -n 1000 -T'])

% Hmm... After normalisation and deconfounding randomise fails:
% --
% Starting permutation 1 (Unpermuted data)
% Warning: The unpermuted statistic image for the current image contains no positive values, and cannot be processed with TFCE. A blank output image will be created.
% --
% Perhaps normalisation violates some of the assumptions of permutation
% testing? Maybe better to model the confounds within randomise itself...

%% Randomise with confounds
OUTDIR = [OUTPUT 'Covid/' FTYPE '_Cov_randomise'];
if ~exist(OUTDIR); mkdir(OUTDIR); end

% Construct the design matrix: [Covid cases, Covid controls, confounds]
X = [CV(AllCovIdx) ~CV(AllCovIdx) conf];

% Copy the contrasts file from the example
tosystem(['cp ' OUTPUT 'Covid/Example_FEAT_setup/cov_vs_con.con ' OUTDIR '/design.con'])

% Write the design.mat file
DESF = [OUTDIR '/design.mat'];
fid = fopen(DESF,'w');
fprintf(fid,'%s\n',['/NumWaves       ' ns(size(X,2))]);
fprintf(fid,'%s\n',['/NumPoints      ' ns(length(AllCovIDs))]);
fprintf(fid,'%s\n','/PPheights              1.000000e+00    1.000000e+00');
fprintf(fid,'%s\n','');
fprintf(fid,'%s\n','/Matrix');
for ii=1:length(AllCovIDs)
    % Find the subject index
    Idx = find(subject_IDs_unique == AllCovIDs(ii));
    
    % Write format: Covid (0=no, 1=yes)   Control (0=no, 1=yes)
    fprintf(fid,'%d\t%d\n',CV(Idx),~CV(Idx));
end
fclose(fid);






%% Later
% Look at Covid/non-Covid, out of interest
% Consider whether we can model out structural effects seen in Gwen's paper
% to make this comparison more valid





%%  ---- Generate figures for ISMRM 2023 abstract ---- %%
% Find an example subject
ID = 1050779;
cd([datdir '2' ns(ID) '/ASL/BASIL/OXASL_ra_dir/']);
!fsleyes qc_output/calib_wholehead.nii.gz qc_output/asldata_diff_mean.nii.gz -dr 0 300 native_space/perfusion_calib.nii.gz -dr 0 100 -cm hot native_space/arrival.nii.gz -dr 0 2.5 -cm hot native_space/aCBV_calib.nii.gz -dr 0 1 -cm blue-green &

% For some reason arrival time map uses a very dilated mask - reapply mask
% here
!fslmaths native_space/arrival.nii.gz -mul native_space/mask.nii.gz /well/okell/projects/biobank/ASL_analysis/tmp_arrival

% Further mask by CBF
!fslmaths native_space/perfusion_calib.nii.gz -thr 25 -bin -mul native_space/arrival.nii.gz /well/okell/projects/biobank/ASL_analysis/tmp_arrival2

%% Example model fit
ASLdat = ra('qc_output/asldata_diff_mean.nii.gz');
CBF = ra('native_space/perfusion_calib');
ATT = ra('native_space/arrival');
M0 = ra('calib/M0');

%FSLvox = [7 28 15]; % Looks too perfect!
%FSLvox = [7 26 15]; 
FSLvox = [9 24 15]; 
Matvox = FSLvox +1;
tmeas = [400:400:2000];
figure; plot(tmeas,squeeze(ASLdat(Matvox(1),Matvox(2),Matvox(3),:)),'o','linewidth',2)
tsim = 300:2100; Alpha = 0.85;
ASLdatsim = M0(Matvox(1),Matvox(2),Matvox(3))*Alpha*BuxtonCASLModel(tsim/1000+1.8,CBF(Matvox(1),Matvox(2),Matvox(3)),ATT(Matvox(1),Matvox(2),Matvox(3)),1.8,1.3,1.65);
hold on;
plot(tsim,ASLdatsim,'-','linewidth',2);
legend('Data','Fit');
xlim([min(tsim) max(tsim)]);
ylim([0 max(ASLdatsim)*1.1])
xlabel 'PLD/s'
ylabel 'ASL signal/au'

toexportfig(gcf,'Example_ASL_fit',FigOutDir,0.2);



%% Look at unique variance explained by ASL IDPs
y = age1(ASLScan1Idx); % Simple age prediction
%y(~ASLScan1Idx) = NaN; % Exclude subjects who don't have ASL data
X = IDPs1_i(ASLScan1Idx,:);
ASLIDPIdx = (IDP_categories==18);

% Univariate associations
[uniR,uniP]=corr(y,X,'rows','pairwise');
IDPNo = 1:length(uniP);
neglogp = -log10(uniP);
scatter(IDPNo(~ASLIDPIdx),neglogp(~ASLIDPIdx));
hold on
scatter(IDPNo(ASLIDPIdx),neglogp(ASLIDPIdx));
xlabel 'IDP number'; ylabel '-log10(p)'; title 'Univariate association with age'
legend('Other IDPs','ASL IDPs')

% Identify the top ~50% 
Top50Idx = neglogp > median(neglogp); % Cut down IDPs to those with biggest univariate associations
NotTooManyNaNsIdx = sum(isnan(X)) < 50; % Remove IDPs with lots of NaNs
Xcut = X(:,Top50Idx & NotTooManyNaNsIdx); % Cut down the design matrix to allow the regression to run
ASLIDPIdxcut = ASLIDPIdx(Top50Idx & NotTooManyNaNsIdx); % Keep track of which remaining IDPs are ASL

% Multi-variate association
[betas,stats]=robustfit(Xcut, nets_normalise(y));

Idx = 1:length(betas); 
figure; 
scatter(Idx(~ASLIDPIdxcut),betas(~ASLIDPIdxcut)); hold on; 
scatter(Idx(ASLIDPIdxcut),betas(ASLIDPIdxcut)); 
legend('Non ASL','ASL')

% From Karla/Steve's NN paper code:
clear poop*; j=1; for i=find(IDP_modality_types(19:end)==3)
  [grot1,grot2]=robustfit(nets_normalise(grotX), nets_normalise(NETdraw(:,i)));
  poopT(j,:)=grot1(2:end);  poopP(j,:)=grot2.p(2:end); j=j+1;
end


% %% Test local functions
% testfn(2)
% loadvars({'age1'},WSPC)
% 
% %% Local functions
% function testfn(blah)
%     disp(blah)
%     testvar = blah;
%     assignin('caller','testvar',testvar)
% end

% Function to load variables from a file into the calling function workspace
% function loadvars(varnames,fname)
%     for ii = 1:length(varnames)
%         disp(['Loading variable: ',varnames{ii}])
%         tmp = load(fname,varnames{ii});
%         disp(['Assigning ',varnames{ii} ' to calling workspace'])
%         assignin('caller',varnames{ii},eval(['tmp.' varnames{ii}]));
%     end
% end