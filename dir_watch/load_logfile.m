%% Load data from a log-file
% 30/01/2017 DK Shin 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Log file format
%-------------------------------
% TDC_RAW1.TXT, DATETIME1, PARAMS11, PARAMS12, ...
% TDC_RAW2.TXT, DATETIME2, PARAMS21, PARAMS22, ...
% ...
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUT
% filedata: a struct with fields: id, datetime, params
%   - id is a positive integer
%   - datetime is a char string
%   - params is a 1D cell-array of strings
%
%-------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function filemeta = load_logfile(path_to_log)
% check file exists
if exist(path_to_log,'file')~=2
    error('File %s does not exist.\nEnter a valid path to log file.',path_to_log);
end

% preallocate output
nlines=linecount(path_to_log);      % num entries in file
filemeta=repmat(struct('id',NaN,'datetime',NaN,'params',NaN),[nlines,1]);

% Read files
fid=fopen(path_to_log,'r');

tline=fgetl(fid);
counter=1;
while ischar(tline)
    telems=strip(strsplit(tline,','));
    telems=telems(1:end-1);     % last elem is empty
    
    % parse line
    tid=dfile2id(telems{1});
    tdatetime=telems{2};        % TODO - currently a char string
    tparams=telems(3:end);
    
    tfilemeta=[];
    tfilemeta.id=tid;
    tfilemeta.datetime=tdatetime;
    tfilemeta.params=tparams;
    
    filemeta(counter)=tfilemeta;
    
    tline=fgetl(fid);
    counter=counter+1;
end
fclose(fid);    % close file

end