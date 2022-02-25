% script for getting the (DRNL) filter parameters (Lopez-Poveda, meddis 2001)
% Author: Morten L?ve Jepsen, 2.nov 2005, rev. 15 feb 2006, 19 feb 2007
%
% usage:  [linDRNLparOut,nlinDRNLparOut] = getDRNLparam(CF);
%
% The returned DRNLparam is a strucures containing the parameter name and
% values

% ------------------------ %
% v 2.0 by Borys K. 23rd April 2015
% v 3.0 by Helia December, 2018
% v 4.0 by Helia March, 2019

function [linDRNLparOut,NlinDRNLparOut] = getDRNLparam_HI(CF,subj);

%% DRNL for normal hearing, Morten 2007

%% Initialize the parameter structures
linDRNLstruct  = struct('parname',{},'vals',{}); % initialize linear paramater vector
NlinDRNLstruct = struct('parname',{},'vals',{}); % initialize nonlinear paramater vector
linDRNLparOut  = struct('parname',{},'vals',{});  NlinDRNLparOut = struct('parname',{},'vals',{});

linDRNLstruct(1).parname = 'CF_lin';  linDRNLstruct(2).parname = 'nGTfilt_lin';
linDRNLstruct(3).parname = 'BW_lin';  linDRNLstruct(4).parname = 'g';
linDRNLstruct(5).parname = 'LP_lin_cutoff';  linDRNLstruct(6).parname = 'nLPfilt_lin';
linDRNLparOut=linDRNLstruct;

NlinDRNLstruct(1).parname = 'CF_nlin';  NlinDRNLstruct(2).parname = 'nGTfilt_nlin';
NlinDRNLstruct(3).parname = 'BW_nlin';  NlinDRNLstruct(4).parname = 'a';
NlinDRNLstruct(5).parname = 'b';  NlinDRNLstruct(6).parname = 'c';
NlinDRNLstruct(7).parname = 'LP_nlin_cutoff';  NlinDRNLstruct(8).parname = 'nLPfilt_nlin';
NlinDRNLparOut = NlinDRNLstruct;

%% Common parameters not subject to immediate change by HI
linDRNLstruct(1).vals = 10^(-0.06762+1.01679*log10(CF)); % Hz, CF_lin,
linDRNLstruct(2).vals = 3; % number of cascaded gammatone filters
linDRNLstruct(3).vals = 10^(.03728+.75*log10(CF)); % Hz, BW_lin.
linDRNLstruct(5).vals = 10^(-0.06762+1.01*log10(CF)); % Hz, LP_lin cutoff
linDRNLstruct(6).vals = 4; % no. of cascaded LP filters,
NlinDRNLstruct(1).vals = 10^(-0.05252+1.01650*log10(CF)); % Hz, CF_nlin
NlinDRNLstruct(2).vals = 3; % number of cascaded gammatone filters,
NlinDRNLstruct(3).vals = 10^(-0.03193+.77*log10(CF)); % Hz, BW_nlin
NlinDRNLstruct(7).vals = 10^(-0.05252+1.01650*log10(CF)); % LP_nlincutoff
NlinDRNLstruct(8).vals = 3; % no. of cascaded LP filters in nlin path,



%%  NH PARAMETERS
if strcmp(subj,'NH') % Model for normal-hearing CASP
    linDRNLstruct(4).vals = 10^(4.20405 -.47909*log10(CF)); %g
    if CF<=1000
        NlinDRNLstruct(4).vals = 10^(1.40298+.81916*log10(CF)); % a,
        NlinDRNLstruct(5).vals = 10^(1.61912-.81867*log10(CF)); % b [(m/s)^(1-c)]
    else
        NlinDRNLstruct(4).vals = 10^(1.40298+.81916*log10(1500)); % a,
        NlinDRNLstruct(5).vals = 10^(1.61912-.81867*log10(1500)); % b [(m/s)^(1-c)]
    end
    NlinDRNLstruct(6).vals = 10^(-.60206); % c, compression coeff
    
    %% Alternative for NH above 1500 Hz
    % if strcmp(subj,'NH') % Model for normal-hearing CASP
    %     linDRNLstruct(4).vals = 10^(4.20405 -.47909*log10(CF)); %g
    %     if CF<=1500 %1000
    %         NlinDRNLstruct(4).vals = 10^(1.40298+.81916*log10(CF)); % a, the 1500 assumption is no good for compressionat low freq filters
    %         NlinDRNLstruct(5).vals = 10^(1.61912-.81867*log10(CF)); % b [(m/s)^(1-c)]
    %     else
    %         NlinDRNLstruct(4).vals = 10^(1.40298+.77*log10(CF)); % a, the 1500 assumption is no good for compressionat low freq filters
    %         NlinDRNLstruct(5).vals = 10^(1.61912-.80*log10(CF)); % b [(m/s)^(1-c)]
    %     end
    %     NlinDRNLstruct(6).vals = 10^(-.60206); % c, compression coeff
    % end
    
    
    %% HI PARAMETERS
elseif strcmp(subj,'HIx') %
    %% HI listener with no compression
    linDRNLstruct(4).vals = 10^(4.20405 -.47909*log10(CF)); %g
    NlinDRNLstruct(4).vals = 0; % a,
    NlinDRNLstruct(5).vals = 10^(1.61912-.81867*log10(CF)); % b [(m/s)^(1-c)]
    NlinDRNLstruct(6).vals = .25; % c, compression coeff
elseif strncmp(subj,'OHC', 3)
    %% Fixed OHCloss profiles -  by Helia R. Iborra 23.1.2017
    
    % Common parameters:
    if CF<=1000
    NlinDRNLstruct(5).vals = 10^(1.61912-.81867*log10(CF)); % b [(m/s)^(1-c)]
    else
     NlinDRNLstruct(5).vals = 10^(1.61912-.81867*log10(1500)); % b [(m/s)^(1-c)]
    end
    
    NlinDRNLstruct(6).vals = .25; % c, compression coeff
    linDRNLstruct(4).vals = 10^(4.20405 -.47909*log10(CF)); %g
    
    if strcmp(subj,'OHCloss_5dB')
        if CF<=1000
            NlinDRNLstruct(4).vals = 10^(1.40298+.81916*log10(CF))*0.562341325; % a,
        else
            NlinDRNLstruct(4).vals = 10^(1.40298+.81916*log10(1500))*0.562341325; % a,
        end
    elseif strcmp(subj,'OHCloss_10dB')
        if CF<=1000
            NlinDRNLstruct(4).vals =10^(1.40298+.81916*log10(CF))*0.316227766; % a,
        else
            NlinDRNLstruct(4).vals = 10^(1.40298+.81916*log10(1500))*0.316227766; % a,
        end
    elseif strcmp(subj,'OHCloss_20dB')
        if CF<=1000
            NlinDRNLstruct(4).vals = 10^(1.40298+.81916*log10(CF))*0.1; % a,
        else
            NlinDRNLstruct(4).vals = 10^(1.40298+.81916*log10(1500))*0.1; % a,
        end
    elseif strcmp(subj,'OHCloss_30dB')
        if CF<=1000
            NlinDRNLstruct(4).vals = 10^(1.40298+.81916*log10(CF))*0.0562341325; % a,
        else
            NlinDRNLstruct(4).vals = 10^(1.40298+.81916*log10(1500))*0.0562341325; % a,
        end
    elseif strcmp(subj,'OHCloss_40dB')
        if CF<=1000
            NlinDRNLstruct(4).vals = 10^(1.40298+.81916*log10(CF))*0.031622777; % a,
        else
            NlinDRNLstruct(4).vals = 10^(1.40298+.81916*log10(1500))*0.031622777; % a,
        end
    end
    
      
        
        
else
 %% Individual HI subjects - read from files (OHC, BMCE)
% Modified by Helia Relano 20.01.2018 - does not require changing the code
% when adding subjects  

[CFdiff,lookupnum] = min(abs([50:25:10000] - CF));
    

% Common parameters:

linDRNLstruct(4).vals = 10^(4.20405 -.47909*log10(CF)); %g

    % Common parameters:
    if CF<=1000
    NlinDRNLstruct(5).vals = 10^(1.61912-.81867*log10(CF)); % b [(m/s)^(1-c)]
    else
     NlinDRNLstruct(5).vals = 10^(1.61912-.81867*log10(1500)); % b [(m/s)^(1-c)]
    end
    
% Individual parameters: 

indivFilename = ['HI-fitting/OHCloss_', subj, '.mat'];
    
if exist(indivFilename) == 2
        load(indivFilename);
        
  if CF<=1000
            NlinDRNLstruct(4).vals = 10^(1.40298+.81916*log10(CF))*out(lookupnum); % a,
     else
            NlinDRNLstruct(4).vals = 10^(1.40298+.81916*log10(1500))*out(lookupnum); % a,
  end
     
    else
        msg = [indivFilename, ' does not exist!'];
        error(msg);
end


    indivFilename = ['HI-fitting/BMCE_', subj, '.mat'];
    
if exist(indivFilename) == 2
        load(indivFilename);
            
     NlinDRNLstruct(6).vals = out(lookupnum); % c, compression coeff
     
else
    
     NlinDRNLstruct(6).vals = 0.25; % If not specified, c = 0.25
     
end    
     
     
% % %     %% Individual HI subjects - read from files
% % %     %% Modified by Borys Kowalewski 01.03.2016 - does not require changing
% % %     % the code when adding subjects
% % %     [CFdiff,lookupnum] = min(abs([50:25:10000] - CF));
% % %     
% % %     %% Check g
% % %     indivFilename = ['models/subjfits/param_g_', subj, '.mat'];
% % %     if exist(indivFilename) == 2
% % %         load(indivFilename);
% % %         linDRNLstruct(4).vals = 10.^(param(lookupnum));
% % %     else
% % %         msg = [indivFilename, ' does not exist!'];
% % %         error(msg);
% % %     end
% % %     
% % %     %% Check a
% % %     indivFilename = ['models/subjfits/param_a_', subj, '.mat'];
% % %     if exist(indivFilename) == 2
% % %         load(indivFilename);
% % %         NlinDRNLstruct(4).vals = 10.^(param(lookupnum));
% % %     else
% % %         msg = [indivFilename, ' does not exist!'];
% % %         error(msg);
% % %     end
% % %     
% % %     %% Check b
% % %     indivFilename = ['models/subjfits/param_b_', subj, '.mat'];
% % %     if exist(indivFilename) == 2
% % %         load(indivFilename);
% % %         NlinDRNLstruct(5).vals = 10.^(param(lookupnum));
% % %     else
% % %         msg = [indivFilename, ' does not exist!'];
% % %         error(msg);
% % %     end
% % %     
% % %     %% Check c
% % %     indivFilename = ['models/subjfits/param_c_', subj, '.mat'];
% % %     if exist(indivFilename) == 2
% % %         load(indivFilename);
% % %         NlinDRNLstruct(6).vals = 10.^(param(lookupnum));
% % %     else
% % %         NlinDRNLstruct(6).vals = 0.25;
% % %         msg = [indivFilename, ' does not exist! c set to 0.25'];
% % %         warning(msg);
% % %     end
% % %     
% % %     %% End of modification by Borys
end

%% Other params
for k=1:6
    linDRNLparOut(k).vals = linDRNLstruct(k).vals;
end
for k=1:8
    NlinDRNLparOut(k).vals = NlinDRNLstruct(k).vals;
end

%% Old fittings by Morten
% if strcmp(subj,'HI1') % fitted 28.01.2009, complete
%     [CFdiff,lookupnum] = min(abs([50:25:10000] - CF));
%     load param_g_HI1.mat; linDRNLstruct(4).vals = 10.^(param(lookupnum));%linDRNLstruct(4).vals = param(lookupnum)
%     load param_a_HI1.mat; NlinDRNLstruct(4).vals = 10.^(param(lookupnum));%NlinDRNLstruct(4).vals = param(lookupnum)
%     load param_b_HI1.mat; NlinDRNLstruct(5).vals = 10.^(param(lookupnum));%NlinDRNLstruct(5).vals = param(lookupnum)
%     NlinDRNLstruct(6).vals = 0.25; % c, compression coeff
% end
% if strcmp(subj,'HI2') % 29.01.2009, complete
%     [CFdiff,lookupnum] = min(abs([50:25:10000] - CF));
%     load param_g_HI2.mat; linDRNLstruct(4).vals = 10.^(param(lookupnum));%linDRNLstruct(4).vals = param(lookupnum)
%     load param_a_HI2.mat; NlinDRNLstruct(4).vals = 10.^(param(lookupnum));%NlinDRNLstruct(4).vals = param(lookupnum)
%     load param_b_HI2.mat; NlinDRNLstruct(5).vals = 10.^(param(lookupnum));%NlinDRNLstruct(5).vals = param(lookupnum)
%     NlinDRNLstruct(6).vals = 0.25; % c, compression coeff
% end
% if strcmp(subj,'HI3') % fitted 29 jan 2009, complete
%     [CFdiff,lookupnum] = min(abs([50:25:10000] - CF));
%     load param_g_HI3.mat; linDRNLstruct(4).vals = 10.^(param(lookupnum));%linDRNLstruct(4).vals = param(lookupnum)
%     load param_a_HI3.mat; NlinDRNLstruct(4).vals = 10.^(param(lookupnum));%NlinDRNLstruct(4).vals = param(lookupnum)
%     load param_b_HI3.mat; NlinDRNLstruct(5).vals = 10.^(param(lookupnum));%NlinDRNLstruct(5).vals = param(lookupnum)
%     load param_c_HI3.mat; NlinDRNLstruct(6).vals = 10.^(param(lookupnum)); % c, compression coeff
% end
% if strcmp(subj,'HI4') % IMK, moderate HL, PTA @1 kHz = 20 dB HL, and @4 kHz = 35 dB HL
%     [CFdiff,lookupnum] = min(abs([50:25:10000] - CF));
%     load param_g_HI4.mat; linDRNLstruct(4).vals = 10.^(param(lookupnum));%linDRNLstruct(4).vals = param(lookupnum)
%     load param_a_HI4.mat; NlinDRNLstruct(4).vals = 10.^(param(lookupnum));%NlinDRNLstruct(4).vals = param(lookupnum)
%     load param_b_HI4.mat; NlinDRNLstruct(5).vals = 10.^(param(lookupnum));%NlinDRNLstruct(5).vals = param(lookupnum)
%     NlinDRNLstruct(6).vals = 0.25; % c, compression coeff
% end
% if strcmp(subj,'HI5') % fitted 30.01 2009, complete
%     [CFdiff,lookupnum] = min(abs([50:25:10000] - CF));
%     load param_g_HI5.mat; linDRNLstruct(4).vals = 10.^(param(lookupnum));%linDRNLstruct(4).vals = param(lookupnum)
%     load param_a_HI5.mat; NlinDRNLstruct(4).vals = 10.^(param(lookupnum));%NlinDRNLstruct(4).vals = param(lookupnum)
%     load param_b_HI5.mat; NlinDRNLstruct(5).vals = 10.^(param(lookupnum));%NlinDRNLstruct(5).vals = param(lookupnum)
%     load param_c_HI5.mat; NlinDRNLstruct(6).vals = 10.^(param(lookupnum)); % c, compression coeff
% end
% if strcmp(subj,'HI6') % fitted 30.01 2009, complete 1 kHz
%     [CFdiff,lookupnum] = min(abs([50:25:10000] - CF));
%     load param_g_HI6.mat; linDRNLstruct(4).vals = 10.^(param(lookupnum));%linDRNLstruct(4).vals = param(lookupnum)
%     load param_a_HI6.mat; NlinDRNLstruct(4).vals = 10.^(param(lookupnum));%NlinDRNLstruct(4).vals = param(lookupnum)
%     load param_b_HI6.mat; NlinDRNLstruct(5).vals = 10.^(param(lookupnum));%NlinDRNLstruct(5).vals = param(lookupnum)
%     NlinDRNLstruct(6).vals = 0.25; % c, compression coeff
% end
% if strcmp(subj,'HI7') % fitted 30.01.2009
%     [CFdiff,lookupnum] = min(abs([50:25:10000] - CF));
%     load param_g_HI7.mat; linDRNLstruct(4).vals = 10.^(param(lookupnum));%linDRNLstruct(4).vals = param(lookupnum)
%     load param_a_HI7.mat; NlinDRNLstruct(4).vals = 10.^(param(lookupnum));%NlinDRNLstruct(4).vals = param(lookupnum)
%     load param_b_HI7.mat; NlinDRNLstruct(5).vals = 10.^(param(lookupnum));%NlinDRNLstruct(5).vals = param(lookupnum)
%     load param_c_HI7.mat; NlinDRNLstruct(6).vals = 10.^(param(lookupnum)); % c, compression coeff
% end
% if strcmp(subj,'HI8') %
%     [CFdiff,lookupnum] = min(abs([50:25:10000] - CF));
%     load param_g_HI8.mat; linDRNLstruct(4).vals = 10.^(param(lookupnum));%linDRNLstruct(4).vals = param(lookupnum)
%     load param_a_HI8.mat; NlinDRNLstruct(4).vals = 10.^(param(lookupnum));%NlinDRNLstruct(4).vals = param(lookupnum)
%     load param_b_HI8.mat; NlinDRNLstruct(5).vals = 10.^(param(lookupnum));%NlinDRNLstruct(5).vals = param(lookupnum)
%     NlinDRNLstruct(6).vals = 0.25; % c, compression coeff
% end
% if strcmp(subj,'HI9') % fitted 30.01.2009
%     [CFdiff,lookupnum] = min(abs([50:25:10000] - CF));
%     load param_g_HI9.mat; linDRNLstruct(4).vals = 10.^(param(lookupnum));%linDRNLstruct(4).vals = param(lookupnum)
%     load param_a_HI9.mat; NlinDRNLstruct(4).vals = 10.^(param(lookupnum));%NlinDRNLstruct(4).vals = param(lookupnum)
%     load param_b_HI9.mat; NlinDRNLstruct(5).vals = 10.^(param(lookupnum));%NlinDRNLstruct(5).vals = param(lookupnum)
%     NlinDRNLstruct(6).vals = 0.25; % c, compression coeff
% end
% if strcmp(subj,'HI10') %
%     [CFdiff,lookupnum] = min(abs([50:25:10000] - CF));
%     load param_g_HI10.mat; linDRNLstruct(4).vals = 10.^(param(lookupnum));%linDRNLstruct(4).vals = param(lookupnum)
%     load param_a_HI10.mat; NlinDRNLstruct(4).vals = 10.^(param(lookupnum));%NlinDRNLstruct(4).vals = param(lookupnum)
%     load param_b_HI10.mat; NlinDRNLstruct(5).vals = 10.^(param(lookupnum));%NlinDRNLstruct(5).vals = param(lookupnum)
%     load param_c_HI10.mat; NlinDRNLstruct(6).vals = 10.^(param(lookupnum)); % c, compression coeff
% end
%% Same fashion, added by Borys
% if strcmp(subj,'tst') %
%     [CFdiff,lookupnum] = min(abs([50:25:10000] - CF));
%     load param_g_tst.mat; linDRNLstruct(4).vals = 10.^(param(lookupnum));%linDRNLstruct(4).vals = param(lookupnum)
%     load param_a_tst.mat; NlinDRNLstruct(4).vals = 10.^(param(lookupnum));%NlinDRNLstruct(4).vals = param(lookupnum)
%     load param_b_tst.mat; NlinDRNLstruct(5).vals = 10.^(param(lookupnum));%NlinDRNLstruct(5).vals = param(lookupnum)
%     load param_c_tst.mat; NlinDRNLstruct(6).vals = 10.^(param(lookupnum)); % c, compression coeff
% end
% if strcmp(subj,'gme') %
%     [CFdiff,lookupnum] = min(abs([50:25:10000] - CF));
%     load param_g_gme.mat; linDRNLstruct(4).vals = 10.^(param(lookupnum));%linDRNLstruct(4).vals = param(lookupnum)
%     load param_a_gme.mat; NlinDRNLstruct(4).vals = 10.^(param(lookupnum));%NlinDRNLstruct(4).vals = param(lookupnum)
%     load param_b_gme.mat; NlinDRNLstruct(5).vals = 10.^(param(lookupnum));%NlinDRNLstruct(5).vals = param(lookupnum)
%     load param_c_gme.mat; NlinDRNLstruct(6).vals = 10.^(param(lookupnum)); % c, compression coeff
% end
% if strcmp(subj,'hpp') %
%     [CFdiff,lookupnum] = min(abs([50:25:10000] - CF));
%     load param_g_hpp.mat; linDRNLstruct(4).vals = 10.^(param(lookupnum));%linDRNLstruct(4).vals = param(lookupnum)
%     load param_a_hpp.mat; NlinDRNLstruct(4).vals = 10.^(param(lookupnum));%NlinDRNLstruct(4).vals = param(lookupnum)
%     load param_b_hpp.mat; NlinDRNLstruct(5).vals = 10.^(param(lookupnum));%NlinDRNLstruct(5).vals = param(lookupnum)
%     NlinDRNLstruct(6).vals = 0.25; % c, compression coeff
% end
% if strcmp(subj,'lgm') %
%     [CFdiff,lookupnum] = min(abs([50:25:10000] - CF));
%     load param_g_lgm.mat; linDRNLstruct(4).vals = 10.^(param(lookupnum));%linDRNLstruct(4).vals = param(lookupnum)
%     load param_a_lgm.mat; NlinDRNLstruct(4).vals = 10.^(param(lookupnum));%NlinDRNLstruct(4).vals = param(lookupnum)
%     load param_b_lgm.mat; NlinDRNLstruct(5).vals = 10.^(param(lookupnum));%NlinDRNLstruct(5).vals = param(lookupnum)
%     NlinDRNLstruct(6).vals = 0.25; % c, compression coeff
% end
% if strcmp(subj,'ajo') %
%     [CFdiff,lookupnum] = min(abs([50:25:10000] - CF));
%     load param_g_ajo.mat; linDRNLstruct(4).vals = 10.^(param(lookupnum));%linDRNLstruct(4).vals = param(lookupnum)
%     load param_a_ajo.mat; NlinDRNLstruct(4).vals = 10.^(param(lookupnum));%NlinDRNLstruct(4).vals = param(lookupnum)
%     load param_b_ajo.mat; NlinDRNLstruct(5).vals = 10.^(param(lookupnum));%NlinDRNLstruct(5).vals = param(lookupnum)
%     NlinDRNLstruct(6).vals = 0.25; % c, compression coeff
% end
% if strcmp(subj,'ilw') %
%     [CFdiff,lookupnum] = min(abs([50:25:10000] - CF));
%     load param_g_ilw.mat; linDRNLstruct(4).vals = 10.^(param(lookupnum));%linDRNLstruct(4).vals = param(lookupnum)
%     load param_a_ilw.mat; NlinDRNLstruct(4).vals = 10.^(param(lookupnum));%NlinDRNLstruct(4).vals = param(lookupnum)
%     load param_b_ilw.mat; NlinDRNLstruct(5).vals = 10.^(param(lookupnum));%NlinDRNLstruct(5).vals = param(lookupnum)
%     NlinDRNLstruct(6).vals = 0.25; % c, compression coeff
% end


%% EOF