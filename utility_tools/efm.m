classdef efm < nltop
    % efm - class for  manipulation of efm signals
    %   Detailed explanation goes here

    properties
        GUID = 'GUID'
        sourceFile= 'fileName'
        createDate=datetime;
        signals = struct ('measure','','sensor','','monitor','monitor', 'segdat',segdat);

    end

    methods
        function e = efm()
            %UNTITLED2 Construct an instance of this class
            %   Detailed explanation goes here

        end
        function EFMout = combineSensors (EFMin)
            % sigOut = emEFMCombineSensors (sigIn) - comnbine signals from different
            % sensors
            %   Combine Hr1 and HR2 into a single Hr2 signale
            versionNum='1.02';
            EFMout=EFMin;
            sigIn=EFMin.signals;
            measureList= {'HR' 'MHR' 'UA'};
            nMeasure=length(measureList);
            kOut=0;
            sigOut=sigIn(1);
            for iMeasure=1:nMeasure,
                jMeasure=find(startsWith({sigIn.measure},measureList{iMeasure}));
                nSeg4Measure=length(jMeasure);
                if nSeg4Measure>0,
                    kOut=kOut+1;
                    sigOut(kOut)=efm.sigCombine(sigIn(jMeasure));
                end
            end
            EFMout.comment= [EFMin.comment '; emEFMCombineSensors ' versionNum];
            EFMout.createDate=date;
            EFMout.signals=sigOut;
        end

        function S=coverage2CTG (e, CTG, reprocessFlag)
            ctgConn=CTG.connection;
            if nargin<3
                reprocessFlag=true
            end
            C=coverage(e);
            [nSamp,nChan,nReal]=size(C);
            cDomain=domain(C);
            cVal=double(C);
            jSamp=0;
            GUID=e.GUID;
            S=struct;
            for iChan=1:nChan

                for iSamp=1:nSamp
                    if cVal(iSamp,iChan)>0
                        jSamp=jSamp+1;
                        S(jSamp).GUID={e.GUID};
                        S(jSamp).measure={e.signals(iChan).measure};
                        S(jSamp).epochStart=cDomain(iSamp);
                        S(jSamp).pctCoverage=cVal(iSamp,iChan);
                    end
                end
            end

            % Check to see if records exist for this GUID
            % only execute if reprocessing
            if reprocessFlag,
                sql=['select * from signalCoverage where GUID= ' '''' GUID ''''];
                T1=fetch(ctgConn,sql);
                % if so delete them
                if ~isempty(T1)
                    sqlDelete =['delete from signalCoverage where GUID= ' '''' GUID ''''];
                    exec(ctgConn, sqlDelete);
                end
            end
            % Write out records for this GUID
            if jSamp>0
                T=struct2table(S);
                sqlwrite(ctgConn,'signalCoverage',T);
            end
        end


        function C = coverage(e, epochLen)
            % Determine coverage of signals in e
            % Epoch length - length of epoch in minutes
            if nargin==1,
                epochLen=20;
            end

            startTime=-12*60*60; % in seconds
            endTime=-0.25;
            domainIncr=0.25;
            nSamp4Epoch=epochLen*60/domainIncr;
            nSec4Epoch=epochLen*60;
            refDomain=[startTime:domainIncr:endTime]';
            nSamp=length(refDomain);
            refSig=-ones(nSamp,1);
            nSig=length(e.signals);
            for iSig=1:nSig,
                curSig=e.signals(iSig).segdat;
                dStart=curSig.domainStart;
                nEpoch=floor(dStart(1)/nSec4Epoch);
                startTime=nEpoch*nSec4Epoch;
                % in seconds
                endTime=-0.25;
                domainIncr=0.25;
                nSamp4Epoch=epochLen*60/domainIncr;
                refDomain=[startTime:domainIncr:endTime]';
                nSamp=length(refDomain);
                refSig=-ones(nSamp,1);
                nSig=length(e.signals);
                curSig=e.signals(iSig).segdat;
                curMeasure=e.signals(iSig).measure;
                chanName{iSig}=curMeasure;
                curNL=nldat(curSig);
                curDomain=domain(curNL);
                curIDX=idx4domain(startTime,domainIncr,curDomain);
                curSig=refSig;
                curSig(curIDX)=double(curNL);
                % find nans

                iBad=find(isnan(curSig));
                curSig(iBad)=-1;
                if strcmp(curMeasure,'UA'),
                    iBad=find(curSig>=127);
                else
                    iBad=find(curSig==0);
                end
                curSig(iBad)=-1;
                iCount=0;
                for i=1:nSamp4Epoch:nSamp,
                    iCount=iCount+1;
                    j=i+nSamp4Epoch-1;
                    k=length(find(curSig(i:j)==-1));
                    C(iCount,iSig)=100*(nSamp4Epoch-k)/nSamp4Epoch;
                end
                C=nldat(C);
                set(C,'chanNames',chanName, 'domainIncr', epochLen, 'domainStart',startTime/60, ...
                    'domainName', 'Epoch start (Minutes)','chanUnits',{' perCent'}, ...
                    'comment',[ e.GUID ': Signal coverage per epoch']);

            end








        end

        function efmStruct= efm2struct ( e, doConvertSegdat )
            if nargin < 2
                doConvertSegdat = false;
            end
            fieldList=fieldnames(e);
            for i=1:length(fieldList)
                curField=fieldList{i};
                efmStruct.(curField)=get(e,curField);
            end
            if doConvertSegdat
                for j=1:length(e.signals)
                    efmStruct.signals(j).segdat=segdat2struct(e.signals(j).segdat);
                end
            end
        end

        function EE=epoch(E, measure,epochStart, epochLen)
            % EE=epoch(E, measure,epochStart, epochLen)
            % Return an epoch of one or more measues in an EFM file
            % measure - cell arry of measures to return
            % epochStart - start time of epoch (minutes)
            % epochLen - length of epoch in (minutes)
            % EE - nldat containing signals
            if epochStart>0,
                error ('epochStart must beless than zero');
            end
           
            if ischar(measure)
                measure={measure};
            end

            epochStart=epochStart*60;
            epochLen=epochLen*60;
            nSig=length(measure);
            s=E.signals;
            domainIncr=E.signals(1).segdat.domainIncr;
            mList={s.measure};
            tMin=-12*60*60;
            tMax=0;
            idx=tMin:domainIncr:tMax;
            epochEnd=epochStart+epochLen-domainIncr;
            epochDomain=epochStart:domainIncr:epochEnd;
            epochIdx=idx4domain(tMin,domainIncr,epochDomain);

            for iSig=1:nSig,
                iMeasure=find(strcmp(measure{iSig},mList));
                if  isempty(iMeasure),
                    error(['Measure not found:' measure{iSig}]);
                end

                seg=E.signals(iMeasure).segdat;

                tempX=nan(length(idx),1);
                segNL=nldat(seg);
                segD=domain(segNL);
                segIdx=idx4domain(tMin,domainIncr,segD);
                tempX(segIdx)=double(segNL);
                epochData(:,iSig)=tempX(epochIdx);
            end
            EE=segNL;
            set(EE,'dataSet',epochData,'domainStart',epochStart,'chanNames',measure, 'comment', ...
                ['epoch from: ' E.GUID]);

        end
        function l = signalEQ(e1,e2)
            % test for equality of sginals in two efm files
            % returns true of data isi the same
            % false if not;
            s1=e1.signals;
            n1=length(s1);
            s2=e2.signals;
            n2=length(s2);
            if n1 ~=n2
                l=false;
                return
            end
            for i=1:n1,
                seg1=s1(i).segdat.dataSet;
                seg2=s2(i).segdat.dataSet;
                if length(seg1) ~= length(seg2)
                    l=false;
                    return
                end
                if ~all(seg1==seg2)
                    l=false;
                    return
                end
                l=true';
            end
        end





        function eSig = getMeasure (eIn, measure )
            % getMeasure (eIn, measure) - returns one or more signal from an efm
            % object
            % eIn - inout efm objvect
            % measure - cell array of one or more measures to return
            %          - measure may be the full name or just the start. An
            %          error occurs more than one meaaure start with the
            %          same name. 
            % Generate an error if a specified measure is not found.
            if ~iscell(measure)
                measure={measure};
            end
            nMeasure=length(measure);
            eSig=eIn;
            sig=eIn.signals;
            measureList={sig.measure};
            ptr=[];
            for iMeasure=1:nMeasure,
                curPtr=find(strcmp(measure{iMeasure},measureList));
                if isempty (curPtr)
                    %error (['Measure not found:' measure{iMeasure}]);
                    continue
                end
                %ptr(iMeasure)=curPtr;
                ptr = [ptr curPtr];
            end
            if isempty(ptr)
                measure
                disp('Measures not found:');
                eSig=[];
                return
            end
            signal=eIn.signals(ptr);
            eSig.signals=signal;
        end
        
        function s = getSegdat(eIn, measure)
            e=getMeasure(eIn, measure);
            s= e.signals.segdat;
        end


        function eInter = intersect ( eIn, mhrFlag)
            % intersect returns interection of signal in signal list
            % Ein - input EFM object
            % Eout- signal to intersect (defaults is HR2 UA)
            % mhrFlag - false, returns intesction of FHR & UA
            % mhrFalg - truye, reutrns intescton of FHR,UA,& MHR
            % eInter - efm object with intsection of signals
            %        - empty f there is no interection
            if nargin==1,
                mhrFlag=false;
            end
            sig=eIn.signals;
            measureList={sig.measure};
            eInter=[];
            s1=getMeasure(eIn,'HR2');
            if isempty(s1)
                disp('No HR signal');
                return
            end

            s2=getMeasure(eIn,'UA');
            if isempty(s2)
                disp('No US signal');
                return
            end
            zSeg=intersect(s1.signals.segdat,s2.signals.segdat);
            if isempty(zSeg)
                eInter=[];
                return
            end

            eInter=efm;
            eInter.signals(1)=s1.signals;
            eInter.signals(1).segdat=segdat(zSeg(:,1));
            eInter.signals(2)=s2.signals;
            eInter.signals(2).segdat=segdat(zSeg(:,2));
            eInter.comment='Intersection: HR2 UA';
            if mhrFlag
                s3=getMeasure(eIn,'MHR');
                if isempty(s3);
                    disp('NO MHR for this trace');
                    return
                end
                zSeg1=intersect(eInter.signals(1).segdat,s3.signals.segdat);
                zSeg2=intersect(eInter.signals(2).segdat,s3.signals.segdat);
                eInter.signals(1)=s1.signals;
                eInter.signals(1).segdat=segdat(zSeg1(:,1));
                eInter.signals(2)=s2.signals;
                eInter.signals(2).segdat=segdat(zSeg2(:,1));
                eInter.signals(3)=s3.signals;
                eInter.signals(3).segdat=segdat(zSeg2(:,2));
                eInter.comment='Intersection of HR2, UA, &MHR';
            end


        end


        function  N=nldat(E)
            % N=nldat(E) convert a EFM object to equivalent nldat objects
            % returns all signals from an EFM object as nldat object with
            % the same start and stop. Missing values are indicated as
            % nans;
            nSig=length(E.signals);
            tMin=0;
            tMax=-12*60*60;
            for iSig=1:nSig,
                d=domain(E.signals(iSig).segdat);
                tMin=min(tMin,min(d));
                tMax=max(tMax,max(d));
            end
            domainIncr=E.signals(1).segdat.domainIncr;
            newDomain=tMin:domainIncr:tMax;
            newData=nan(length(newDomain),nSig);
            chanNameList={};
            for iSig=1:nSig,
                curNL=nldat(E.signals(iSig).segdat);
                curIdx=idx4domain(tMin, domainIncr,domain(curNL));
                newData(curIdx,iSig)=double(curNL);
                chanNameList{iSig}=E.signals(iSig).measure;
            end
            N=nldat(newData,'domainStart',tMin,'domainIncr',domainIncr,'chanNames',chanNameList);
            set(N,'comment', ['nldat of:' E.GUID]);
        end


        function plot (EFM)
            % emEFMPlot - plot a EFdata set
            clf;
            nChan=length(EFM.signals);
            % Determine a comomon start end end time
            for iChan=1:nChan,
                d=domain(EFM.signals(iChan).segdat);
                startTime(iChan)=min(d)/3600;
                endTime(iChan)=max(d)/3600;
            end
            startTime=min(startTime);
            endTime=max(endTime);


            for iChan=1:nChan,
                curEFM=EFM.signals(iChan).segdat;
                set(curEFM, 'domainIncr',curEFM.domainIncr/3600, 'domainName','Hours', ...
                    'domainStart', curEFM.domainStart/3600);
                subplot (nChan,1,iChan);
                plot(curEFM);
                titleStr=strrep([EFM.signals(iChan).measure ' ' EFM.signals(iChan).sensor  ...
                    ' ' EFM.signals(iChan).monitor],'_','\_');
                title(titleStr);

                set(gca,'xlim',[startTime endTime]);
                ylabel(EFM.signals(iChan).measure)
                xlabel('');
            end
            xlabel('Hours');
            streamer(EFM.GUID);
        end


        function sigOut = repair (sigIn, nBadMax, nInterp)
            %% sigOut = repair (sigIn, nBadMax, nInterp)
            % repair an efm object by interpolating missing/bad values.
            %  Long Periods of bad/missing data result in creation of segment.
            %
            % nBadmMax [50]  - maxmimum number of sequential bad samples to interpolate.
            % nInterp [5]   - Number of points before and after bad segmenet to use for interpolation
            if nargin <3
                nBadMax=50;  % Maximum length of a bad sample event to interpolate
            end
            if nargin <2
                nInterp=5;
            end
            sigOut=sigIn;
            set(sigOut,'comment','efm.repair');
            nSig=length(sigIn.signals);
            for iSig=1:nSig

                curMeasure=sigIn.signals(iSig).measure;
                curSig=sigIn.signals(iSig).segdat;
                curSeg=nldat(curSig);
                newSeg=curSeg;
                curSegLen=length(curSeg);
                c=categorical;
                c(1:curSegLen)='good';
                iNan=find(isnan(curSeg));
                c(iNan)='bad' ;
                switch curMeasure
                    case {'HR1' 'HR2' 'MHR'}
                        iZero=find(double(curSeg)==0);
                        c(iZero)='bad';
                    case 'UA'
                        iZero=find(double(curSeg)>=127);
                        c(iZero)='bad';
                    otherwise
                        disp(['No code for measure:' curMeasure]);
                end
                e=eseq.cseq2eseq(c);
                % Find bad events

                % Convert 2 bad events separated by a short good event to one long bad
                % event
                iBad=find([e.type]=='bad');
                nBad=length(iBad);
                for jBad=1:nBad
                    kBad=iBad(jBad);
                    iNext=kBad+1;
                    if iNext>length(e)
                        break
                    end
                    if e(iNext).type == 'good' & e(iNext).nSamp <nInterp
                        e(iNext).type='bad';
                    end
                end
                c=cseq(e);
                e=eseq.cseq2eseq(c);
                iBad=find([e.type]=='bad');
                nBad=length(iBad);

                %% interpolate values for short bad events
                for jBad=1:nBad   % Loop for bad events in seq
                    kBad=iBad(jBad);
                    iBstart=e(kBad).startIdx;
                    iBend=e(kBad).endIdx;
                    iFront=iBstart-(nInterp:-1:1);
                    iEnd=iBend+[1:nInterp];
                    xSpline=cat(2,iFront,iEnd);
                    xInterp=iBstart:iBend;
                    if e(kBad).nSamp<nBadMax % short bad event so interpolate

                        % Handle specal cases
                        %  Bad event is at end of record so set to nans
                        if any(xSpline>curSegLen) | any(xSpline<1)
                            newSeg(xInterp)=nan;
                        else
                            ySpline=double(curSeg(xSpline))';
                            yInterp=round(spline(xSpline,ySpline,xInterp));
                            newSeg(xInterp)=yInterp';
                        end
                    else % Long event so set values to nans
                        newSeg(xInterp)=nan;
                    end
                end
                newSegDat=segdat(newSeg);
                sigOut.signals(iSig).segdat=newSegDat;
            end
        end







        function  save (EFM,  EFMtype, CTG, processComment )
            % emEFMsave (EFM, saveDir, EFMtype, conn)
            % save EFM to saveDir. If file already exists it is overwritten
            %
            GUID=EFM.GUID;
            saveLocation=efm.getLocation(EFMtype, GUID);
            i=strfind(saveLocation,'\');
            dirName=saveLocation(1:(max(i)-1));
            if ~exist(dirName)
                disp(['Creating: ' dirName]);
                cmd1=['!mkdir ' dirName];
                eval(cmd1);
            end

            %saveLocation='X:\output\';
            ctgConn=CTG.connection;


            disp(['Saving: ' saveLocation]);
            efmStruct=efm2struct(EFM);
            save(GUID, 'efmStruct');
            cmd=['!move /Y ' GUID '.mat ' saveLocation ];
            eval(cmd);

            % Update database with results
            % Check to see if record(s) exists for this phase, if so delete it
            sql=['select * from CTGstatus where processPhase=' '''' EFMtype '''' ...
                ' AND GUID=' '''' GUID '''' ];
            Texist=fetch(ctgConn,sql);
            for i=1:height(Texist),
                sql=['delete from CTGstatus where recordId =' num2str(Texist.recordId(i))];
                execute(ctgConn,sql)
            end
            %% Insert new record
            T=table;
            T.GUID={GUID};
            T.sourceFile={EFM.sourceFile};
            T.processPhase={EFMtype};
            T.processVersion={EFM.comment};
            T.processDate=datetime;
            T.processComment={processComment};
            sqlwrite(ctgConn,'CTGstatus', T);
        end
        function SS = stats(EFM)
            % chanStats = emEFMstats(EFM) - compute important statistics of a EM EFM file
            % input:
            %   XS - a struct generated by emReadFile
            % output:
            %  chanStats - a structure array continaing statistics of each segment for al segdat objects in XS
            % fields include:
            %  segNum - segment number
            %  segLen - length of segment
            %  startTime - time segment starts (in hours before birth)
            %  endTime- time at which segment ends(in hours before birth)
            %  nDropout - number of dropouts (HR=0 or UA >- 127
            %  pctDropout - percentage of segment samples that are dropouts
            %  nDropEvent - number of dropout Events
            %  dropOutLength
            %  Quantiles - 50,75 and 9% quantiles of dropoutl lenghts


            %% Segment statistics
            chanStats=struct;
            SS=struct;
            SS.GUID={EFM.GUID};
            SS.comment={'emEFMstats V1.01'};
            SS.createDate=datetime;
            iSeg=0;
            signal=EFM.signals;
            nChan=length(signal);
            for iChan=1:nChan,
                XS=signal(iChan).segdat;
                nSeg=segCount(XS);
                curMeasure={signal(iChan).measure};
                SS.signal(iChan).measure={signal(iChan).measure};
                SS.signal(iChan).sensor={signal(iChan).sensor};
                SS.signal(iChan).monitor={signal(iChan).monitor};
                SS.signal(iChan).nSeg=nSeg;
                d=domain(XS);
                SS.signal(iChan).startTime= min(d)/3600;
                SS.signal(iChan).endTime=max(d)/3600;
                SS.signal(iChan).nSamp=length(XS);
                SS.signal(iChan).pctCoverage=100*length(XS)/length(nldat(XS));
                curData=double(signal(iChan).segdat);
                [ nBad, pctBad, nBadEvent] =efm.badDataStats (curMeasure, curData);
                SS.signal(iChan).nBad=nBad;
                SS.signal(iChan).pctBad=pctBad;
                SS.signal(iChan).nBadEvent=nBadEvent;
                %


                for iSeg=1:nSeg,

                    SS.signal(iChan).segment(iSeg).segNum=iSeg;
                    curSeg=segGet(XS,iSeg);
                    curData=double(curSeg);
                    nSamp=length(curSeg);
                    SS.signal(iChan).segment(iSeg).nSamp=nSamp;
                    d=domain(curSeg);
                    SS.signal(iChan).segment(iSeg).startTime= min(d)/3600;
                    SS.signal(iChan).segment(iSeg).endTime=max(d)/3600;
                    %         signal(iChan.segment(iSeg).numDuplicateSampleTimes=SS.signal(iChan).numDuplicateSampleTimes;
                    %         signal(iChan.segment(iSeg).nullValCnt=SS.signal(iChan).nullValCnt;
                    %          signal(iChan.segment(iSeg).eqValCnt=SS.signal(iChan).eqValCnt;
                    %         signal(iChan.segment(iSeg).diffValCnt=SS.signal(iChan).diffValCnt;
                    [ nBad, pctBad, nBadEvent] =efm.badDataStats (curMeasure, curData);
                    SS.signal(iChan).segment(iSeg).nBad=nBad;
                    SS.signal(iChan).segment(iSeg).pctBad=pctBad;
                    SS.signal(iChan).segment(iSeg).nBadEvent=nBadEvent;
                    %

                end
            end

        end

        function disp(sys)
            builtin('disp',sys);
            s=sys.signals;
            nSig=length(s);
            for iSig=1:nSig
                disp(['Signal ' num2str(iSig) ':' s(iSig).measure]);
            end
        end



    end

    methods(Static)

        function e = loadGUID (GUID, EFMtype, CTG)
            %  e = loadGUID (GUID, EFMtype) load an EFM file for a GUID
            % GUID - guid of file to load
            % EFMtype - type of EFM file ('combined'/'repair');
            %

            if nargin==1,
                EFMtype='combined';
            end
            if strcmp(EFMtype,'repair'),
                S=patternsStatus(CTG,GUID);
                if isempty(S)
                    i=[];
                else
                    i=find( strcmp( 'EFM from PatternsMerge',S.status_str));
                end
                if isempty(i)
                    disp(['GUID not found: ' GUID]);
                    e=[];
                    return
                else
                    fileLocation=S.dest_file{i};
                end
            else
                fileLocation= efm.getLocation(EFMtype,GUID);
            end
            if ~exist(fileLocation),
                disp(['File not found: ' fileLocation]);
                e=[];
            else
                e=efm.loadFile(fileLocation);
            end

        end


        function e = loadFile (fileName)
            % e = loadFile (fileName)
            % load an efm file
            load(fileName)
            if exist('efmStruct')
                e=efm.struct2efm(efmStruct);
            elseif exist('EFM')
                if isstruct(EFM),
                    e=efm.struct2efm(EFM);
                elseif isa(EFM,'efm')
                    e=EFM;
                else
                    error (['Bad format in file:' fileName]);
                end
            end
        end

        function multiPlot (g)
            % multiPlot (g) - loads and plots multiple EFM files in the same window
            % g = cell array of guids to compare
            % currently only works for two guids
            e1=efm.loadGUID(g{1});
            e2=efm.loadGUID(g{2});
            % compare multile EFM files
            measureList={'HR2' 'UA'};
            for i=1:2,
                curMeasure=measureList{i};
                m1=nldat(getMeasure(e1,curMeasure));
                m2=nldat(getMeasure(e2,curMeasure));
                subplot (2,1,i);
                if length(m1.dataSet)>0
                    h=line(m1); set(h,'color','b');
                end
                if length(m2.dataSet)>0
                    h=line(m2);
                    set(h,'color','r');
                end
                title(curMeasure);
                if i==1,
                    legend(g);
                end
            end
        end

        function  EFM=readFile (fileName, sensorsToIgnore)
            % function  [ZS]=emEFMreadFile (fileName)
            % Read in a compressed EFM csv file, parse it, return values
            % Input:
            %   fileName - full file specification for file
            % Output:
            % ZS - a structure array with one entry for each unique measure/sensor/monitor combination.
            % ZS has fields:
            %   comment
            % measure
            %   sensor
            %   monitor
            %   segdat - parsed as a nlid_tools/segdat object
            %   other fields - a bunch of other fields to test for multiple samples.
            %%
            if nargin==1,
                sensors2Ignore= {'INOP' 'No_Trans' 'No-Trans'}; % List of sensors to ignore
            end
            versionNum='2.01]';
            f=unzip(fileName,pwd);
            i=max(strfind(fileName,'\'));
            GUID=fileName(i+1:end);
            j=strfind(GUID,'.');
            GUID=GUID(1:j-1);
            % GUID=strrep(GUID,'.zip','');
            %[numVal,strVal,allVal]=xlsread(f{1});
            inputT=readtable(f{1});
            nRow=height(inputT);
            %% Initialization
            if nargin==1,
                sensorsToIgnore={};
            end
            domainStart=[];
            onsetPointer=[];
            segLength=0;
            duplicateSampleFlag=false;
            emptySampleFlag=false;
            nullValCnt=0;        %Count of duplicate time samples with one null and one good value
            eqValCnt=0;         % Count of dubplicate time samples with the same values;
            diffValCnt=0;       %Count of duplicate time samples with different values;
            doubleNullCnt=0;   % Count of dubpcate time samples where both are null
            %
            %% Column definitions
            timeInCol=1;
            timeRecCol=2;
            measurementCol=3;
            monitorCol=4;
            sensorCol=5;
            readingCol=6;

            %% parse channel information
            X=struct('T',[]);
            chanNameCol=2;
            sampleCol=5 ; % determine full channel name
            % check to see if there is more than one monitor
            uniqueMonitor=unique(inputT.monitor(:));
            nMonitor=length(uniqueMonitor);
            % strVal=strrep(strVal,'-','_');  % change "-" to "_" so it can be sued as a channel name
            % allChanNames={};
            % for iRow=2:nRow
            %     curChanName=cat (2, strVal{iRow,sensorCol}, '_', strVal{iRow,measurementCol});
            %     if ~ismember(curChanName, allChanNames),
            %         allChanNames=cat(1,allChanNames, curChanName);
            %     end
            % end
            % chanNameList=unique(allChanNames);
            % nChan=length(chanNameList);
            %
            % Initialize output structure
            measure=inputT.measurement(:);
            monitor=inputT.monitor(:);
            sensor=inputT.sensor(:);
            timeVal=(inputT.deltaInMillisecs(:)/1000.); % Convert to milliseconds
            timeRecVal=(inputT.deltaRecordedInMillisecs(:)/1000.); % Convert to milliseconds
            sampVal=(inputT.reading(:));
            uniqueMeasure=unique(measure);
            nMeasure=length(uniqueMeasure);
            %% Parse data
            ZS={};
            iVar=0;
            for iMonitor=1:nMonitor,
                curMonitor=uniqueMonitor{iMonitor};
                for iMeasure=1:nMeasure,
                    curMeasure=uniqueMeasure{iMeasure};
                    iCurMeasure=find(strcmp(curMonitor,monitor) & strcmp(curMeasure, measure));
                    sensor4measure=unique(sensor(iCurMeasure));
                    nSensor=length(sensor4measure);
                    for iSensor=1:nSensor,
                        curSensor=sensor4measure{iSensor};
                        curChanName=cat(2,curMeasure, '_', curSensor);
                        curChanName=cat(2,curChanName,'_',curMonitor);
                        is4m = find(strcmp(curMonitor,monitor) & strcmp(curMeasure, measure) & strcmp(curSensor,sensor));
                        curSampVal=sampVal(is4m,:);
                        curTimeVal=timeVal(is4m);
                        duplicateSampleTimes=[];
                        emptySampleTimes=[];
                        nullValCnt=0;        %Count of duplicate time samples with one null and one good value
                        eqValCnt=0;         % Count of dubplicate time samples with the same values;
                        diffValCnt=0;       %Count of duplicate time samples with different values;
                        doubleNullCnt=0 ;  % Count of dubpcate time samples where both are null
                        %

                        %% Handle Duplicate Times if they exist
                        dt=diff(curTimeVal);
                        iDuplicate=find(dt==0);
                        nDuplicate=length(iDuplicate);

                        if nDuplicate>0,
                            duplicateSampleTimes=curTimeVal(iDuplicate);;
                            % disp([num2str(nDuplicate) ' duplicate time samples exist']);
                            iDrop=[];
                            for jDuplicate=1:nDuplicate,
                                iCur=iDuplicate(jDuplicate);
                                iNext=iCur+1;
                                curSampleValid=~isnan(curSampVal(iCur)) & ~isempty(curSampVal(iCur));
                                nextSampleValid=  ~isnan(curSampVal(iNext)) & ~isempty(curSampVal(iNext));

                                % Cuurent  value is a valid entry - drop the next one
                                if curSampleValid & ~nextSampleValid
                                    iDrop=cat(1,iDrop,iNext);
                                    nullValCnt=nullValCnt+1;
                                    % currentsample is not valid but the next is - drop current value
                                elseif ~curSampleValid  & nextSampleValid
                                    iDrop=cat(1,iDrop,iCur);
                                    nullValCnt=nullValCnt+1;
                                    % Both samples valid
                                elseif curSampleValid & nextSampleValid
                                    if curSampVal(iCur)==curSampVal(iNext),
                                        eqValCnt=eqValCnt+1;
                                    else
                                        diffValCnt=diffValCnt+1;
                                    end
                                    iDrop=cat(1,iDrop,iNext);
                                else % Neither sample valid
                                    doubleNullCnt=doubleNullCnt+1;
                                    iDrop=cat(1,iDrop,iCur,iNext);
                                end
                            end
                            curSampVal(iDrop)=[];
                            curTimeVal(iDrop)=[];
                        end

                        %% Get rid of remaining samples with nans
                        iNan=find(isnan(curSampVal));
                        nEmpty=length(iNan);

                        if length(iNan)>0,
                            emptySampleTimes=curTimeVal(iNan);
                            % disp([ num2str(length(iNan)) ' empty samples exist']);
                            curSampVal(iNan)=[];
                            curTimeVal(iNan)=[];
                        end

                        if length(curSampVal)>0,
                            if ismember(curSensor,sensorsToIgnore),
                                disp(['Ignoring sensor:' curSensor]);

                            else

                                Z=nldat(curSampVal,'domainValues',-curTimeVal,'chanNames', ...
                                    {curChanName},'domainIncr',.25,'comment',[GUID ':' curMonitor]);

                                ZStmp=segdat(Z);
                                nSeg=segCount(ZStmp);
                                segInfo={}; % Mis-spelled in V2.0 resultinng in extraneous entries in come cases
                                for i=1:nSeg;
                                    segInfo{i}=[curSensor '-' curMonitor];
                                end
                                ZStmp.segInfo=segInfo;
                                iVar=iVar+1;
                                ZS.GUID=GUID;
                                ZS.sourceFile=fileName;
                                ZS.comment=['Creaed by emReadFile ' versionNum];
                                ZS.createDate=date;
                                ZS.signals(iVar).measure=curMeasure;
                                ZS.signals(iVar).sensor=curSensor;
                                ZS.signals(iVar).monitor=curMonitor;
                                %                     ZS(iVar).emptySampleTimes=emptySampleTimes;
                                %                     ZS(iVar).numEmptySampleTimes=length(emptySampleTimes);
                                %                     ZS(iVar).duplicateSampleTimes=duplicateSampleTimes;
                                %                     ZS(iVar).numDuplicateSampleTimes=length(duplicateSampleTimes);
                                %                     ZS(iVar).nullValCnt=nullValCnt;
                                %                     ZS(iVar).eqValCnt=eqValCnt;
                                %                     ZS(iVar).diffValCnt=diffValCnt;
                                %                     ZS(iVar).doubleNullCnt=doubleNullCnt;
                                ZS.signals(iVar).segdat=ZStmp;
                            end
                        end
                    end
                end
            end
            EFM=efm.struct2efm(ZS);
        end
        function sOut = sigCombine (sIn)
            % combine signals of the same variable from different
            % sensor
            nSig=length(sIn);
            if nSig==1,
                sOut=sIn;
                return
            end

            % Determine which signals are high quality
            % sort so high quaility monitors come last
            sensorList={sIn(:).sensor};
            [sensorListSorted, iSort]=sort(sensorList);
            iSort=flip(iSort);
            sortedSensors=sensorList(iSort);

            %% Sepcial handling when combined H1 ansd H2
            for i=1:nSig
                if contains(sIn(i).measure,'HR')
                    sIn(i).sensor= [ sIn(i).measure '_' sIn(i).sensor]; % HR1 or HR2 to sensor
                    % add H1 or HR2 to segIfo
                    curSegDat=sIn(i).segdat;
                    segInfo=curSegDat.segInfo;
                    nSeg=length(segInfo);
                    for j=1:nSeg
                        segInfo{j}=strrep([ sIn(i).measure '_' segInfo{j}],' ','');
                    end
                    sIn(i).segdat.segInfo=segInfo;
                end
            end
            sSort=sIn(iSort);
            sOut=sSort(1);

            % Concatonate signals
            sCat=sSort(1).segdat;

            for iSig=2:nSig,
                sCat=cat(1,sCat,sSort(iSig).segdat);
            end
            sensorList=[sSort(1).measure '_' sSort(1).sensor];
            monitorList=sSort(1).monitor;
            for i=2:nSig,
                sensorList=[ sensorList '&' sSort(i).measure '_' sSort(i).sensor];
                monitorList=[ monitorList '&' sSort(i).monitor];
            end
            sOut.sensor=sensorList;
            sOut.monitor=monitorList;
            sOut.segdat=sCat;
        end


        function [ nBad, pctBad, nBadEvent] =badDataStats (curMeasure, curData)
            % baddata  stats
            localMeasure=curMeasure{1:2};
            if strcmp('HR', localMeasure) | strcmp('MH', localMeasure),
                xBad=(double(curData)==0); % true if a drop out
            elseif strcmp('UA', localMeasure),
                xBad=(curData>=127);
            else
                xBad=[];
            end
            nBad=length(find(xBad));
            if nBad==0,
                nBadSamp=nBad;
                nBadEvent=0;
                pctBad=0;
            else
                xCat=categorical(xBad);
                xEvent=eseq(xCat);
                nBadEvent=length(find([xEvent.type]=='true'));
                nBadSamp=nBad;
                nBadEvent=nBadEvent;
                pctBad=100*nBad/length(curData);
            end
        end

        function e= struct2efm ( efmStruct )
            e=efm;
            if isempty(efmStruct)
                e.comment='No data';
            else
                fieldList={'GUID' 'sourceFile', 'comment' 'createDate' 'signals'};
                for i=1:length(fieldList)
                    curField=fieldList{i};
                    set(e,curField,efmStruct.(curField));
                end
            end
        end



        function [d,g,f]=efmDir (folder, subFolderFlag)
            % [d,g]=efmDir (folder)
            % Returns directory and GUID information for efm files in
            % folder and subfiler
            % d - directory information
            % g - cell array of GUIDs
            % f - full fileame
            % folder
            if nargin==0,
                folder=efm.getLocation('combined');
            end
            if nargin==1,
                subFolderFlag=false;
            end
            if subFolderFlag
                curDir=pwd;
                cd(folder);
                d=dir('**/*.mat');
                cd(curDir);
            else
                d=dir(folder);
            end
            d=d(3:end);
            n={d.name};
            g=strrep(n,'.mat','');
            nFile=length(d);
            for i=1:nFile,
                f{i}=[ folder '\' d(i).name];
            end
            disp(['Folder: ' folder ' has  ' num2str(length(n)) ' mat files']);
        end

        function d = getLocation(option, spec )
            % efm.getLocaion (optin,spec) returns directory location
            % options are:
            % getLocation('base') - returnss fully qualified pointer to the base dir
            % getLocation('combined') - base directory for combined dEFM files
            % getLocation('raw', N) - directory for raw data set N
            % getLocaton('repair') direcotry for patterns repaired data set
            %
            %
            baseDir= 'X:\';
            switch option
                case 'base'
                    d=baseDir;
                case 'dataset'
                    dataSetTxt=num2str(spec);
                    if length(dataSetTxt)<2
                        dataSetTxt=[ '0' dataSetTxt];
                    end
                    %  d= [ efm.getLocation('raw') 'dataset' dataSetTxt '\' ];
                    d=[ baseDir 'aprilresultssmb\' 'dataset' dataSetTxt '\' ];
                case 'combined'
                    if nargin==1
                        d=[ baseDir 'processing\combined\'];
                    elseif nargin==2
                        guid=char(spec);
                        subDir=guid(1:2);
                        d=[ baseDir 'processing\combined\' subDir '\' guid '.mat'];
                    end
                case 'raw'
                    d=[baseDir 'aprilresults\'];
                case 'rerun'
                    d=[baseDir 'processing\rerun']
                otherwise
                    error (['Option not defined:' option])
            end
        end
        function readAudit(fileList)
            % readAudit(fileList) - compares results fo readFile to raw
            % trace.
            % Generates one figures for each signal that  superimposes
            % the segdat signal (in color) generated by readfile and raw trace
            % data in black.
            if ischar(fileList),
                fileList={fileList};
            end
            nFile=length(fileList);
            % Validate eading against raw data
            for iFile= 1:nFile,
                fileName=fileList{iFile};


                %% Get Raw Data

                R = efm.readFile(fileName);

                %% Get trace dsata and read in

                i=max(strfind(fileName,'\'));
                GUID=fileName(i+1:end);
                j=strfind(GUID,'.');
                GUID=GUID(1:j-1);
                T=readtable([GUID '.csv']);

                %% Compare signals
                nSig=length(R.signals);
                for iSig=1:nSig
                    curSig=R.signals(iSig);
                    iCur=find (strcmp(T.measurement, curSig.measure) ...
                        & strcmp(T.sensor,curSig.sensor) ...
                        & strcmp(T.monitor,curSig.monitor));
                    curReading=T.reading(iCur);
                    %curT=-T.deltaRecordedInMillisecs(iCur)/(1000);
                    curT=-T.deltaInMillisecs(iCur)/(1000);
                    figure(iSig);clf
                    plot(curSig.segdat);

                    h=line(curT,curReading); set(h,'color','k');
                    title([ GUID ' ' curSig.measure ' ' curSig.sensor ' ' curSig.monitor]);
                end
                continueFlag=input_l('Continue',true);
                if ~continueFlag
                    break
                end
            end
        end
    end
end
