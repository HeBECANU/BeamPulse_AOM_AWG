function filemeta = load_logfile(path_to_log)
%Load data from a log-file
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Log file format
%-------------------------------
% TDC_RAW1.TXT, DATETIME1, PARAMS11, PARAMS12, ...
% TDC_RAW2.TXT, DATETIME2, PARAMS21, PARAMS22, ...
% ...
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUT
%   filedata: a 1D struct-array with fields: id, datetime, params
%       - id is a positive integer
%       - datetime is a char string
%       - params is a 1D cell-array of strings
%
%-------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 30/01/2017
% DK Shin 

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

%% Bad attempts from log
% % pad missing shots with repeated param (log file reports only the SUCCESSFUL
% % shot)
% fileid=vertcat(filemeta(:).id);
% id_start=fileid(1);
% id_end=fileid(end);
% 
% nshot=id_end-id_start+1;
% filemeta_padded=repmat(struct('id',NaN,'datetime',NaN,'params',NaN),[nshot,1]);
% counter=nshot;      % count down
% tfilemeta=struct('id',NaN,'datetime',NaN,'params',NaN);
% for ii=id_end:-1:id_start
%     idxInLog=find(fileid==ii);
%     if ~isempty(idxInLog)
%         tfilemeta=filemeta(idxInLog);   % trying this param set
%         filemeta_padded(counter)=tfilemeta;
%     else
%         % pad it with previously set param set
%         filemeta_padded(counter).id=counter;
%         filemeta_padded(counter).params=tfilemeta.params;
%     end
%     counter=counter-1;
% end
% filemeta=filemeta_padded;

end