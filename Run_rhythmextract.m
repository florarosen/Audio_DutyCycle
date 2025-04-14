%% Get file names to input to the RHYTHM_EXTRACT_FFT function 
% Should load all postprocessed files 
clear

genpath = '/Users/florarosenzweig/Dropbox/Mon Mac (mac-114-318.local)/Documents/MATLAB/PROJECTS/XP_Audio_DutyCycles/output/YOUNG/EEG/';
PostprocPath = fullfile(genpath, '/Postprocessing');

%add letswave
rnb_tools_path = '/Users/florarosenzweig/Dropbox/Mon Mac (mac-114-318.local)/Documents/MATLAB/Toolboxes/rnb_tools/src'; 
lw_path = '/Users/florarosenzweig/Dropbox/Mon Mac (mac-114-318.local)/Documents/MATLAB/Toolboxes/letswave6-master'; 


addpath(genpath(rnb_tools_path)); 
addpath(genpath(lw_path)); 

%get list of all processed files rerefed to mastoids
cd(PostprocPath)
list_fft_M1M2 = dir('**/scrpt*fft*reref*M1M2*sub0*EEG*.mat');
list_fft_AVG  = dir('**/scrpt*fft*reref*AVG*sub0*EEG*.mat');

%could add same for reref avg
for i=1:size(list_fft_M1M2,1)
    filename = list_fft_M1M2(i).name;
    disp(filename)
    RHYTHM_EXTRACT_FFT(filename)
end

% %could add same for reref avg
% for i=1:size(list_fft_AVG,1)
%     filename = list_fft_AVG(i).name
%     RHYTHM_EXTRACT_FFT(filename)
% end