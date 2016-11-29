
% analyze data from trigger timestamped data from biopac
% 07/07/2015 - Namanh Kapur - created
% 07/08/2015 - Namanh Kapur - edited
% 07/14/2015 - RRS - added HRV bar graph
% 07/15/2015 - RRS - added pupillometry analysis
% 07/24/2015 - Namanh Kapur - removed NaNs and interpolated gaps
% 07/26/2015 - Namanh Kapur - added bounded box analysis
% 08/04/2015 - Namanh Kapur - added group analysis
% 09/01/2015 - RRS - added analysis for skin temperature

% define params
maxBPM = 250; %maximum beats per minute
before = 2; %seconds
after = 5;  %seconds
freq = 200; %Hz
baseline = 1; %second
numTrials = 20;

% channel mapping
% 02Sat = 1;
GSRChannel = 2;
skinTempChannel = 3;
PulseChannel = 4;
HRChannel = 5;
% HRStatus = 6;
CocaineChannel = 7;
NaturalChannel = 8;
ScrambledChannel = 9;

% analyze pupillometry data
pupilFreq = 1000;
% pupilFreq = 250;

% read in all Biopac Data Files
fileName = 'BiopacData/*.mat';
files = dir(fullfile(pwd,fileName));
numSubjects = size(files,1);

% read in all Eye Tracking Files
pupilName = 'ETData/*.txt'; 
pupilFiles = dir(fullfile(pwd,pupilName));

% initialize data arrays for biopac
tNaturalGSR = zeros(numTrials, before*freq + after*freq, numSubjects);
tCocaineGSR = zeros(numTrials, before*freq + after*freq, numSubjects);

tNaturalTemp = zeros(numTrials, before*freq + after*freq, numSubjects);
tCocaineTemp = zeros(numTrials, before*freq + after*freq, numSubjects);

tNaturalHR = zeros(numTrials, before*freq + after*freq, numSubjects);
tCocaineHR = zeros(numTrials, before*freq + after*freq, numSubjects);

tNaturalPeakTimeDiffs = cell(numTrials, 1, numSubjects);
tCocainePeakTimeDiffs = cell(numTrials, 1, numSubjects);

% initialize arrays
tNaturalPupil = zeros(numTrials, before*pupilFreq + after*pupilFreq, numSubjects);
tCocainePupil = zeros(numTrials, before*pupilFreq + after*pupilFreq, numSubjects);
tScrambledCocainePupil = zeros(numTrials, before*pupilFreq + after*pupilFreq, numSubjects);
tScrambledNaturalPupil = zeros(numTrials, before*pupilFreq + after*pupilFreq, numSubjects);

% create time vector
times = -before:1/freq:after; 
times = times(1:end-1);

% start analysis 
for nSubjects = 1:size(files,1)
    
    biopacFile = ['BiopacData/', files(nSubjects).name]; 
    biopacData = load(biopacFile);
    pupilFile = ['ETData/', pupilFiles(nSubjects).name];
    pupil = readtable(pupilFile, 'Delimiter', 'tab', 'TreatAsEmpty', '.');
    
    %% pupillometry variables
    % find when images were displayed

    allImagesInds = strfind(pupil.SAMPLE_MESSAGE, '.jpg');
    allImagesInds = find(not(cellfun('isempty', allImagesInds)));

    psIndsCocaine = strfind(pupil.SAMPLE_MESSAGE, 'ps_cocaine');
    psIndsCocaine = find(not(cellfun('isempty', psIndsCocaine)));

    psIndsNatural = strfind(pupil.SAMPLE_MESSAGE, 'ps_natural');
    psIndsNatural = find(not(cellfun('isempty', psIndsNatural)));

    naturalInds = strfind(pupil.SAMPLE_MESSAGE, 'natural_');
    naturalInds = find(not(cellfun('isempty', naturalInds)));
    naturalInds = setdiff(naturalInds,psIndsNatural);

    cocaineInds = strfind(pupil.SAMPLE_MESSAGE, 'cocaine_');
    cocaineInds = find(not(cellfun('isempty', cocaineInds)));
    cocaineInds = setdiff(cocaineInds,psIndsCocaine);

    % get all the ps inds
    psInds = strfind(pupil.SAMPLE_MESSAGE, 'ps_');
    psInds = find(not(cellfun('isempty', psInds)));

    psNaturalOrder = zeros(length(psIndsNatural), 1);
    for x = 1:length(naturalInds)
        psNaturalOrder(x) = find(psInds == psIndsNatural(x));
    end

    psCocaineOrder = zeros(length(psIndsCocaine), 1);
    for x = 1:length(naturalInds)
        psCocaineOrder(x) = find(psInds == psIndsCocaine(x));
    end
    
    % data vars
    GSR = biopacData.data(:,GSRChannel);
    HR = biopacData.data(:,HRChannel);
    PULSE = biopacData.data(:,PulseChannel);
    TEMP = biopacData.data(:,skinTempChannel);

    % detrend GSR
    meanGSR = mean(GSR); 
    GSR = detrend(GSR) + meanGSR;
    
    % detrend TEMP
    meanTEMP = mean(TEMP);
    TEMP = detrend(TEMP) + meanTEMP;
    
    % find NN peak times of pulse 
    [~, pulsePeakInds] = findpeaks(PULSE, 'MinPeakHeight', 0, 'MinPeakDistance', 1/(maxBPM/60*freq));
    
    % find times when images were presented
    tNatural = find(diff(biopacData.data(:, NaturalChannel)) > 0);
    tCocaine = find(diff(biopacData.data(:, CocaineChannel)) > 0);
    tScrambledAll = find(diff(biopacData.data(:, ScrambledChannel)) > 0);
    tScrambledNatural = tScrambledAll(psNaturalOrder);
    tScrambledCocaine = tScrambledAll(psCocaineOrder);
    
    %% analyze to throw out trials

    analyzePupillometry = 1; 
    NotFixating = find(isnan(pupil.RIGHT_INTEREST_AREA_ID));

    if (length(NotFixating)/length(pupil.RIGHT_INTEREST_AREA_ID)) > 0.3
        analyzePupillometry = 0;
    end
    
    if analyzePupillometry == 1
        
        % throw out trials that don't follow fixation bounds 80% of the time
        for x = 1:length(naturalInds)
            %finds before and after sample indeces
             sampleBefore = naturalInds(x) - before*pupilFreq;
             sampleAfter = naturalInds(x) + after*pupilFreq;
            %find NaNs during analysis indeces as defined above
             tempPupil = pupil.RIGHT_INTEREST_AREA_ID(sampleBefore:sampleAfter-1); 
             tempNotFixating = find(isnan(tempPupil));
            %choose whether to throw trial
            if length(tempNotFixating) /length(tempPupil) > 0.2
                naturalInds(x) = NaN;
                tNatural(x) = NaN;
            end
        end

        for x = 1:length(cocaineInds)
            %finds before and after sample indeces
             sampleBefore = cocaineInds(x) - before*pupilFreq;
             sampleAfter = cocaineInds(x) + after*pupilFreq;
            %find NaNs during analysis indeces as defined above
             tempPupil = pupil.RIGHT_INTEREST_AREA_ID(sampleBefore:sampleAfter-1); 
             tempNotFixating = find(isnan(tempPupil));
            %choose whether to throw trial
            if length(tempNotFixating) /length(tempPupil) > 0.2
               cocaineInds(x) = NaN;
               tCocaine(x) = NaN;
            end
        end

        for x = 1:length(psIndsCocaine)
            %finds before and after sample indeces
             sampleBefore = psIndsCocaine(x) - before*pupilFreq;
             sampleAfter = psIndsCocaine(x) + after*pupilFreq;
            %find NaNs during analysis indeces as defined above
             tempPupil = pupil.RIGHT_INTEREST_AREA_ID(sampleBefore:sampleAfter-1); 
             tempNotFixating = find(isnan(tempPupil));
            %choose whether to throw trial
            if length(tempNotFixating) /length(tempPupil) > 0.2
                psIndsCocaine(x) = NaN;
                tScrambledCocaine(x) = NaN;
                tCocaine(x) = NaN;
            end
        end

        for x = 1:length(psIndsNatural)
            %finds before and after sample indeces
             sampleBefore = psIndsNatural(x) - before*pupilFreq;
             sampleAfter = psIndsNatural(x) + after*pupilFreq;
            %find NaNs during analysis indeces as defined above
             tempPupil = pupil.RIGHT_INTEREST_AREA_ID(sampleBefore:sampleAfter-1); 
             tempNotFixating = find(isnan(tempPupil));
            %choose whether to throw trial
            if length(tempNotFixating) /length(tempPupil) > 0.2
                psIndsNatural(x) = NaN;
                tScrambledNatural(x) = NaN;
                tNatural(x) = NaN;
            end
        end   
    end
    
    %% interpolate between blinks
    nBlinks = 0;
    nDeleteBlink = 50;
    nSamp = 100;

    blinkOnsetInds = find(diff(pupil.RIGHT_IN_BLINK) == 1);
    blinkOffsetInds = find(diff(pupil.RIGHT_IN_BLINK) == -1);

    if length(blinkOnsetInds) == length(blinkOffsetInds)
        nBlinks = length(blinkOnsetInds);
    elseif length(blinkOnsetInds) > length(blinkOffsetInds) % throw away last blink if at end
        pupil = pupil(1:blinkOnsetInds(end)-1,:);
        blinkOnsetInds(end) = [];
        nBlinks = length(blinkOnsetInds);
    end

    for n = 1:nBlinks

        % remove data before and after blink
        pupil.RIGHT_PUPIL_SIZE(blinkOnsetInds(n)-nDeleteBlink:blinkOnsetInds(n)) = NaN; 

        beforeBlinkInds = blinkOnsetInds(n)+1-nSamp:blinkOnsetInds(n)-nDeleteBlink-1;
        afterBlinkInds = blinkOffsetInds(n)+1+nDeleteBlink:blinkOffsetInds(n)+nSamp;

        beforeBlinks = pupil.RIGHT_PUPIL_SIZE(beforeBlinkInds);
        afterBlinks = pupil.RIGHT_PUPIL_SIZE(afterBlinkInds);
        tVal = [blinkOnsetInds(n)+1-nSamp:blinkOffsetInds(n)+nSamp];

        t = [beforeBlinkInds afterBlinkInds];
        data = [beforeBlinks; afterBlinks]'; 

        p = polyfit(t, data,1);
        pEval = polyval(p, tVal); % interpolated pupilSize values

    %       p = spline(t,data);
    %       pEval = ppval(p, tVal);
          pupil.RIGHT_PUPIL_SIZE(blinkOnsetInds(n)-nDeleteBlink:blinkOffsetInds(n)+nDeleteBlink) = pEval(nSamp-nDeleteBlink:end-nSamp+nDeleteBlink);
    end
    
    %% interpolate any remaining nans now
    if any(isnan(pupil.RIGHT_PUPIL_SIZE))
        nBlinks = 0;
        nDeleteBlink = 50;
        nSamp = 100;

        pupilNans = find(isnan(pupil.RIGHT_PUPIL_SIZE));
        pupilNans = [0; pupilNans];

        blinkInds = find(diff(pupilNans) > 1);
        blinkOnsetInds = pupilNans(blinkInds + 1);
        blinkOffsetInds = pupilNans(blinkInds(2:end) - 1);
        blinkOffsetInds = [blinkOffsetInds; pupilNans(end)];

        if length(blinkOnsetInds) == length(blinkOffsetInds)
            nBlinks = length(blinkOnsetInds);
        elseif length(blinkOnsetInds) > length(blinkOffsetInds) % throw away last blink if at end
            pupil = pupil(1:blinkOnsetInds(end)-1,:);
            blinkOnsetInds(end) = [];
            nBlinks = length(blinkOnsetInds);
        end

        for n = 1:nBlinks

            % remove data before and after blink
            pupil.RIGHT_PUPIL_SIZE(blinkOnsetInds(n)-nDeleteBlink:blinkOnsetInds(n)) = NaN; 

            beforeBlinkInds = blinkOnsetInds(n)+1-nSamp:blinkOnsetInds(n)-nDeleteBlink-1;
            afterBlinkInds = blinkOffsetInds(n)+1+nDeleteBlink:blinkOffsetInds(n)+nSamp;

            beforeBlinks = pupil.RIGHT_PUPIL_SIZE(beforeBlinkInds);
            afterBlinks = pupil.RIGHT_PUPIL_SIZE(afterBlinkInds);
            tVal = [blinkOnsetInds(n)+1-nSamp:blinkOffsetInds(n)+nSamp];

            t = [beforeBlinkInds afterBlinkInds];
            data = [beforeBlinks; afterBlinks]'; 

            p = polyfit(t, data,1);
            pEval = polyval(p, tVal); % interpolated pupilSize values

        %       p = spline(t,data);
        %       pEval = ppval(p, tVal);
              pupil.RIGHT_PUPIL_SIZE(blinkOnsetInds(n)-nDeleteBlink:blinkOffsetInds(n)+nDeleteBlink) = pEval(nSamp-nDeleteBlink:end-nSamp+nDeleteBlink);
        end    
    end
    
    %% detrend pupillometry
    meanPupil = mean(pupil.RIGHT_PUPIL_SIZE); 
    pupil.RIGHT_PUPIL_SIZE = detrend(pupil.RIGHT_PUPIL_SIZE) + meanPupil;

    % create time vector
    pupilTimes = -before:1/pupilFreq:after; 
    pupilTimes = pupilTimes(1:end-1);

    %% analyze trials
    if analyzePupillometry == 1

        % natural
        for x = 1:length(naturalInds)
            if isnan(naturalInds(x))
                tNaturalPupil(x,:,nSubjects) = NaN;
            else
                %finds before and after sample indeces
                 sampleBefore = naturalInds(x) - before*pupilFreq;
                 sampleAfter = naturalInds(x) + after*pupilFreq;
                %finds original pupil size Values before psc analysis
                 tempPupil = pupil.RIGHT_PUPIL_SIZE(sampleBefore:sampleAfter-1);
                %computes mean baseline
                 meanBaselinePupil = mean(pupil.RIGHT_PUPIL_SIZE(sampleBefore:sampleBefore+baseline*pupilFreq));
                %calculates percent signal change
                 calcpscPupil = (((tempPupil/meanBaselinePupil)*100)-100); 
                %sets final GSR array to wanted psc analysis output
                 tNaturalPupil(x,:,nSubjects) = calcpscPupil;
            end
        end

        % cocaine
        for x = 1:length(cocaineInds)
            if isnan(cocaineInds(x))
                tCocainePupil(x,:,nSubjects) = NaN;
            else
                %finds before and after sample indeces
                 sampleBefore = cocaineInds(x) - before*pupilFreq;
                 sampleAfter = cocaineInds(x) + after*pupilFreq;
                %finds original pupil size Values before psc analysis
                 tempPupil = pupil.RIGHT_PUPIL_SIZE(sampleBefore:sampleAfter-1);
                %computes mean baseline
                 meanBaselinePupil = mean(pupil.RIGHT_PUPIL_SIZE(sampleBefore:sampleBefore+baseline*pupilFreq));
                %calculates percent signal change
                 calcpscPupil = (((tempPupil/meanBaselinePupil)*100)-100); 
                %sets final GSR array to wanted psc analysis output
                 tCocainePupil(x,:,nSubjects) = calcpscPupil;
            end
        end

        % scrambled cocaine
        for x = 1:length(psIndsCocaine)
            if isnan(psIndsCocaine(x))
                tScrambledCocainePupil(x,:,nSubjects) = NaN;
            else
                %finds before and after sample indeces
                 sampleBefore = psIndsCocaine(x) - before*pupilFreq;
                 sampleAfter = psIndsCocaine(x) + after*pupilFreq;
                %finds original pupil size Values before psc analysis
                 tempPupil = pupil.RIGHT_PUPIL_SIZE(sampleBefore:sampleAfter-1);
                %computes mean baseline
                 meanBaselinePupil = mean(pupil.RIGHT_PUPIL_SIZE(sampleBefore:sampleBefore+baseline*pupilFreq));
                %calculates percent signal change
                 calcpscPupil = (((tempPupil/meanBaselinePupil)*100)-100); 
                %sets final GSR array to wanted psc analysis output
                 tScrambledCocainePupil(x,:,nSubjects) = calcpscPupil;
            end
        end

        % scrambled natural
        for x = 1:length(psIndsNatural)
            if isnan(psIndsNatural(x))
                tScrambledNaturalPupil(x,:,nSubjects) = NaN;
            else
                %finds before and after sample indeces
                 sampleBefore = psIndsNatural(x) - before*pupilFreq;
                 sampleAfter = psIndsNatural(x) + after*pupilFreq;
                %finds original pupil size Values before psc analysis
                 tempPupil = pupil.RIGHT_PUPIL_SIZE(sampleBefore:sampleAfter-1);
                %computes mean baseline
                 meanBaselinePupil = mean(pupil.RIGHT_PUPIL_SIZE(sampleBefore:sampleBefore+baseline*pupilFreq));
                %calculates percent signal change
                 calcpscPupil = (((tempPupil/meanBaselinePupil)*100)-100); 
                %sets final GSR array to wanted psc analysis output
                 tScrambledNaturalPupil(x,:,nSubjects) = calcpscPupil;
            end
        end

        % subtract off psc of ps images
        for x = 1:length(psIndsCocaine)
            tCocainePupil(x,:,nSubjects) = tCocainePupil(x,:,nSubjects) - tScrambledCocainePupil(x,:,nSubjects);
        end

        % subtract off psc of ps images
        for x = 1:length(psIndsNatural)
            tNaturalPupil(x,:,nSubjects) = tNaturalPupil(x,:,nSubjects) - tScrambledNaturalPupil(x,:,nSubjects);
        end
    end
    
    %% analyze biopac data and convert to PSC compared to baseline
    % natural
    for x = 1:length(tNatural)
        if isnan(tNatural(x))
            tNaturalGSR(x,:,nSubjects) = NaN;
            tNaturalTemp(x,:,nSubjects) = NaN;
            tNaturalHR(x,:,nSubjects) = NaN;
            tNaturalPeakTimeDiffs{x,:,nSubjects} = NaN;
        else
            %finds before and after sample indeces
             sampleBefore = tNatural(x) - before*freq;
             sampleAfter = tNatural(x) + after*freq;
            %grab data
             tempGSR =  GSR(sampleBefore:sampleAfter-1);
             tempTEMP = TEMP(sampleBefore:sampleAfter-1);
             tempHR = HR(sampleBefore:sampleAfter-1);
            %compute mean baseline
             meanBaselineGSR = mean(GSR(sampleBefore:sampleBefore+baseline*freq));
             meanBaselineTEMP = mean(TEMP(sampleBefore:sampleBefore+baseline*freq));
             meanBaselineHR = mean(HR(sampleBefore:sampleBefore+baseline*freq));
            %calculate percent signal change
             calcpscGSR = (((tempGSR/meanBaselineGSR)*100)-100); 
             calcpscTEMP = (((tempTEMP/meanBaselineTEMP)*100)-100); 
             calcpscHR = (((tempHR/meanBaselineHR)*100)-100); 
            %sets final GSR array to wanted psc analysis output
             tNaturalGSR(x,:,nSubjects) = calcpscGSR;
             tNaturalTemp(x,:,nSubjects) = calcpscTEMP;
             tNaturalHR(x,:,nSubjects) = calcpscHR;
            %grab peak time diffs for HRV
             trialPeakInds = find(pulsePeakInds >= tNatural(x) & pulsePeakInds <= sampleAfter);
             peakDiffs = diff(pulsePeakInds(trialPeakInds))/freq;
             tNaturalPeakTimeDiffs{x,:,nSubjects} = peakDiffs; 
        end
    end

    % cocaine
    for x = 1:length(tCocaine)
        if isnan(tCocaine(x))
            tCocaineGSR(x,:,nSubjects) = NaN;
            tCocaineTemp(x,:,nSubjects) = NaN;
            tCocaineHR(x,:,nSubjects) = NaN;
            tCocainePeakTimeDiffs{x,:,nSubjects} = NaN;
        else        
            %finds before and after sample indeces
             sampleBefore = tCocaine(x) - before*freq;
             sampleAfter = tCocaine(x) + after*freq;
            %finds original GSR and HR Values before psc analysis
             tempGSR = GSR(sampleBefore:sampleAfter-1);
             tempTEMP = TEMP(sampleBefore:sampleAfter-1);
             tempHR = HR(sampleBefore:sampleAfter-1);
            %computes mean baseline
             meanBaselineGSR = mean(GSR(sampleBefore:sampleBefore+baseline*freq));
             meanBaselineTEMP = mean(TEMP(sampleBefore:sampleBefore+baseline*freq));
             meanBaselineHR = mean(HR(sampleBefore:sampleBefore+baseline*freq));
            %calculates percent signal change
             calcpscGSR = (((tempGSR/meanBaselineGSR)*100)-100); 
             calcpscTEMP = (((tempTEMP/meanBaselineTEMP)*100)-100); 
             calcpscHR = (((tempHR/meanBaselineHR)*100)-100); 
            %sets final GSR array to wanted psc analysis output
             tCocaineGSR(x,:,nSubjects) = calcpscGSR;
             tCocaineTemp(x,:,nSubjects) = calcpscTEMP;
             tCocaineHR(x,:,nSubjects) = calcpscHR;
            %grab peak time diffs for HRV
             trialPeakInds = find(pulsePeakInds >= tCocaine(x) & pulsePeakInds <= sampleAfter);
             peakDiffs = diff(pulsePeakInds(trialPeakInds))/freq;
             tCocainePeakTimeDiffs{x,:,nSubjects} = peakDiffs;
        end
    end
    
    
end

%% plot PSTH 
plotGroupPSTH(tNaturalGSR, tCocaineGSR, freq, 2, 'Natural GSR', 'Cocaine GSR', -1, 1);
plotGroupPSTH(tNaturalTemp, tCocaineTemp, freq, 2, 'Natural Temp', 'Cocaine Temp', -0.1, 0.1);
plotGroupPSTH(tNaturalHR, tCocaineHR, freq, 2, 'Natural HR', 'Cocaine HR', -1, 1);
plotGroupPSTH(tNaturalPupil, tCocainePupil, pupilFreq, 2, 'Natural Pupil', 'Cocaine Pupil', -15, 15);

%% plot bar graphs for HRV
subjNaturalHRV = zeros(size(tNaturalPeakTimeDiffs,3),1);
subjCocaineHRV = zeros(size(tCocainePeakTimeDiffs,3),1);

for x = 1:length(subjNaturalHRV)
    subjNaturalHRV(x) = nanstd(cell2mat(tNaturalPeakTimeDiffs(:,:,x)));
    subjCocaineHRV(x) = nanstd(cell2mat(tCocainePeakTimeDiffs(:,:,x)));
end

figure; 
bar(1, mean(subjNaturalHRV), 'facecolor', 'g');
hold on;
bar(2, mean(subjCocaineHRV), 'facecolor', 'r');

errorbar(1, mean(subjNaturalHRV), nanstd(subjNaturalHRV)/sqrt(length(subjNaturalHRV)))
errorbar(2, mean(subjCocaineHRV),nanstd(subjCocaineHRV)/sqrt(length(subjCocaineHRV)))

scatter(linspace(0.75,1.25,length(subjNaturalHRV)), subjNaturalHRV, 'ko');
scatter(linspace(1.75,2.25,length(subjCocaineHRV)), subjCocaineHRV, 'ko');

set(gca, 'XTick', [1 2]);
set(gca, 'XTickLabel', {'Natural', 'Cocaine'});
set(gca,'YTick', [0 0.25 0.5])
ylim([0 0.5])
set(gca,'FontSize',20);
ylabel('HRV [sec]', 'FontSize', 20);