function D = flbio( fid, option, D )
% FLBIO - perform io operations on flb files
%
%   option:
%       'read_details - read details of next case
%       'read_data' - read data for next case
%       'read_case'
%       'read_index'
%       'read_all'
%       'write_case'
%       'new_real'
%
% Usage
%
% Example
%
% See also
%
%% AUTHOR: Robert Kearney
%% $DATE: 02-Feb-2006 08:48:32 $
%% $Revision: 1.4 $
%% DEVELOPED: 7.1.0.246 (R14) Service Pack 3
%% FILENAME: flbio.m

% DL made changes to lines:75,87,88

switch lower(option),
    case 'read_details'
        D=flbrdet (fid);
    case 'read_case'
        D=flbrdet (fid);
        if isempty(D),
            if feof(fid),
                error('End of file');
            else
                error('Error reading file');
            end
        end
        D=flbrdata(fid, D);
    case 'read_index'
        caseNum=0;
        D=struct([]);
        while 1,
            tempD =flbrdet(fid);
            if isempty(tempD),
                break
            else
                nByte = tempD.nChan * tempD.chanLen *tempD.nReal * tempD.chanFormat;
                status = fseek(fid, nByte, 'cof');
                caseNum=caseNum+1;
                if caseNum==1,
                    D=tempD;
                else
                    D(caseNum)=tempD;
                end
            end
        end
    case 'read_all'
        caseNum=0;
        while 1,
            tempD =flbrdet(fid);
            if isempty(tempD),
                break
            else
                tempD=flbrdata(fid, tempD);
                caseNum=caseNum+1;
                D(caseNum)=tempD;
            end
        end
        
    case 'skip_case'
        D=flbSkipCase(fid);
    case 'write_case'
        dataSize=size(D.Data);
        % DL changed nDim -> D.nDim
        D.nDim=length(dataSize);
        
        % Number of realizations
        if D.nDim==3,
            D.nReal=dataSize(3);
        else
            D.nReal=1;
        end
        D.nChan=dataSize(2);
        D.chanLen=dataSize(1);
        %DL changed chanMax/Min -> D.chanMax/Min
        D.chanMax=max(D.Data);
        D.chanMin=min(D.Data);
        s=flbwdet(fid,D);
        count = fwrite(fid, D.Data, 'float');
        D=1;
    case 'new_real'
        caseIdx=flbio(fid,'read_index');
        nCase=length(caseIdx);
        frewind(fid);
        for i=1:nCase-1,
            s=flbio(fid,'skip_case');
        end
        casePos=ftell(fid);
        oldDetails= flbio(fid,'read_details');
        % Check that details match
        [nSamp, nChan, nReal]=size(D.Data);
        if nSamp ~= oldDetails.chanLen,
            error ('Realization has a different sample length');
        end
        if nChan ~= oldDetails.nChan,
            error ('Realization has a differnet number channels');
        end
        %% Update Details
        Details=oldDetails;
        Details.nDim=3;
        Details.nReal=Details.nReal+ nReal;
        fseek(fid, casePos, 'bof');
        s=flbwdet(fid, Details);
        fseek(fid,0,'eof');
        count = fwrite(fid, D.Data, 'float');
        D=1;
        
    otherwise
        error (['flbio option not defined: ' option ]);
end

function D = flbrdet( fid)
% FLBRDET ...
%  Read next set of details from a flb file
%   ...
% Usage
%
% Example
%
% See also
% read details from a flb file
flb_version = fread(fid,1,'int32');
if isempty(flb_version),
    D=[];
    return
end

if (flb_version < 2),
    disp ('Header error')
    num_chan=-2;
    return
end
if flb_version== 2,
    D.nDim=2;
    D.nReal=1;
    D.nChan=fread(fid,1,'int32')
    D.chanLen=fread(fid,1,'int32')
    D.chanFormat=fread(fid,1,'int32')
    strLen=fread(fid,1,'int32')
    D.domainName=char(fread(fid,strLen,'char'))'
    D.domainIncrement=fread(fid,1,'real*4')
    D.domainStart=fread(fid,1,'real*4')
    commentLen=fread(fid,1,'int32')
    D.comment=char(fread(fid,commentLen,'char'))'
    for iChan=1:D.nChan,
        nameLen=fread(fid,1,'int32')
        D.chanName{iChan}=char(fread(fid,nameLen,'char'))'
    end
    D.chanMin=fread(fid,D.nChan,'real*8')
    D.chanMax=fread(fid,D.nChan,'real*8')
else
    D.nDim=fread(fid,1,'int32');
    D.nReal=fread(fid,1,'int32');
    D.nChan=fread(fid,1,'int32');
    D.chanLen=fread(fid,1,'int32');
    D.chanFormat=fread(fid,1,'int32');
    stringLen=fread(fid,1,'int32');
    D.domainName=char(fread(fid,stringLen,'char'))';
    D.domainIncr = fread(fid,1,'single');
    D.domainStart = fread(fid,1,'single');
    stringLen=fread(fid,1,'int32');
    D.comment=char(fread(fid,stringLen,'char'))';
    
    for iChan=1:D.nChan,
        stringLen=fread(fid,1,'int32');
        D.chanName{iChan}=char(fread(fid,stringLen,'char'))';
    end
    D.chanMin =fread(fid,D.nChan,'float64');
    D.chanMax =fread (fid,D.nChan,'float64');
end
return

%%

function D= flbrdata (fid, D)
% Read next data set
nSamp = D.nChan * D.chanLen *D.nReal;

switch D.chanFormat
    case 2
        formatSpec='short';
    case 4
        formatSpec='float';
end
tempC= fread(fid, nSamp, formatSpec);
C = reshape (tempC, D.chanLen, D.nChan, D.nReal);
D.Data=C;
return

%% Skip case

function status = flbSkipCase (fid)
D=flbrdet(fid);
nByte = D.nChan * D.chanLen *D.nReal * D.chanFormat;
status = fseek(fid, nByte, 'cof');
return

%% Write details

function status = flbwdet( fid, D)
% FLBRDET ...
%  Write a case to current file location
%  D is a structre with fields:
%       D.comment
%       D.chanName
%       D.domainIncr
%       D.domainName
%       D.domainStat
%       D.Datahelp
%   ...
% Usage
%
% Example
%
% See also
% read details from a flb file
flbVersion =  4;
count= fwrite(fid,flbVersion,'int32');
%% Number of dimensions

count = fwrite(fid,D.nDim,'int32');
count=fwrite(fid,D.nReal,'int32');
% number of channels
count=fwrite(fid,D.nChan,'int32');
% Channel Length
count=fwrite(fid,D.chanLen,'int32');
% Data Format
count=fwrite(fid,4,'int32');
% DomainName
stringLen=length(D.domainName);
count=fwrite(fid, stringLen, 'int32');
count=fwrite(fid, D.domainName,'char');
% Domain Increment
count=fwrite (fid,D.domainIncr,'single');
count=fwrite (fid,D.domainStart,'single');
% Comment
stringLen=length(D.comment);
count=fwrite(fid,stringLen,'int32');
count=fwrite(fid,D.comment,'char');
% Channel Names
for iChan=1:D.nChan,
    stringLen=length(D.chanName{iChan});
    count=fwrite(fid,stringLen,'int32');
    count=fwrite(fid, D.chanName{iChan},'char');
end

count=fwrite(fid,D.chanMin,'float64');
count=fwrite(fid,D.chanMax,'float64');


%% data



status=1;
return













