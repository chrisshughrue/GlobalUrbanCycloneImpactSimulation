%%%This code simulates losses to production in urban
%%%areas and subsequent economic cascades through forward and backward
%%%industrial shortages
%%%Author: Chris Shughrue
%%%Oct 15, 2019

%%%Please cite (Shughrue at al, 2020) paper when using this code.

%%%Run parameters
outputPath = '\\simulationRuns\\'
inputPath = '\\inputs\\'
trackInd = 1;%index of historic storm track to simulate. This cycles over all storm indices (1:1209)
paramNum = 0;%switch for sensitivity analysis run
paramVal = 1;%parameter for sensitivity analysis

%load cyclone impact data
load(sprintf('%scycloneImpacts.mat',inputPath))
stormTrackInds = find(sum(wsMat,1)>0);
trackInd = stormTrackInds(trackInd);

%identify existing simulation output
dirl = dir(inputPath);
fileCondition = isempty(find(strcmp({dirl.name},sprintf('runOutputStorm%dParameter%dTest%d.mat',trackInd,paramNum,paramVal))));

if fileCondition==1
    save(sprintf('%srunOutputStorm%dParameter%dTest%d.mat',outputPath,trackInd,paramNum,paramVal),'fileCondition');
    load(sprintf('%senvironmentVariables.mat',inputPath))
    xSumInds = ones(numSectors,1)*[1:numCities]; xSumInds = xSumInds(:);
    sRun = uStar;
    ttMax = 365*5;

    incomingDemandAnomaly = zeros(numel(cityVectIndTrunc),ttMax);
    outgoingDemandAnomaly = zeros(numel(cityVectIndTrunc),ttMax);
    incomingDemandAnomalyCorr = zeros(numel(cityVectIndTrunc),ttMax);
    outgoingDemandAnomalyCorr = zeros(numel(cityVectIndTrunc),ttMax);

    %adjust variables based on input
    if paramNum==1
        ALPHAMAXSENSITIVITY = paramVal;
        TAUSENSITIVITY = 1;
        ELASTICITYSENSITIVITY = 1;
        GAMMASENSITIVITY = 1;
        REBUILDSENSITIVITY = 1;
    elseif paramNum==2
        ALPHAMAXSENSITIVITY = 1;
        TAUSENSITIVITY = paramVal;
        ELASTICITYSENSITIVITY = 1;
        GAMMASENSITIVITY = 1;
        REBUILDSENSITIVITY = 1;
    elseif paramNum==3
        ALPHAMAXSENSITIVITY = 1;
        TAUSENSITIVITY = 1;
        ELASTICITYSENSITIVITY = paramVal;
        GAMMASENSITIVITY = 1;
        REBUILDSENSITIVITY = 1;
    elseif paramNum==4
        ALPHAMAXSENSITIVITY = 1;
        TAUSENSITIVITY = 1;
        ELASTICITYSENSITIVITY = 1;
        GAMMASENSITIVITY = paramVal;
        REBUILDSENSITIVITY = 1;
    elseif paramNum==5
        ALPHAMAXSENSITIVITY = 1;
        TAUSENSITIVITY = 1;
        ELASTICITYSENSITIVITY = 1;
        GAMMASENSITIVITY = 1;
        REBUILDSENSITIVITY = paramVal;
    else
        ALPHAMAXSENSITIVITY = 1;
        TAUSENSITIVITY = 1;
        ELASTICITYSENSITIVITY = 1;
        GAMMASENSITIVITY = 1;
        REBUILDSENSITIVITY = 1;
    end

    %initialize some variables
    jD = jDStar;
    jDGamma = ones(size(jD));
    alphaMax = 1.2.*ALPHAMAXSENSITIVITY;
    overallTau = 180.*TAUSENSITIVITY;
    weeklyTau = 180.*TAUSENSITIVITY;
    elasticity = 0.9.*ELASTICITYSENSITIVITY;
    demandElasticity = 0.9.*ELASTICITYSENSITIVITY;
    gammaP = (0.07./30)*GAMMASENSITIVITY;
    rebuildTau = 1300*REBUILDSENSITIVITY;

    %reset some variables
    tCount = 1;
    T = tStar;
    lambdaX = zeros(numCities*numSectors,1);
    profitDeltaMin = 0.001;
    profitDelta = 0.001 + (rand(size(lambdaX)).*0.001);
    profitFrac = ones(size(lambdaX));
    profitFracTMinus = 0.5 + rand(size(profitFrac));
    profitCalcTMinus = 0.5 + rand(size(profitFrac));
    totalXVect = zeros(numel(cityVectIndTrunc),ttMax);
    totalProfitVect = zeros(numel(cityVectIndTrunc),ttMax);
    limitationMat = zeros(numel(cityVectIndTrunc),ttMax,4);
    limInds = repmat(xSumInds,1,4) + ((ones(numSectors*numCities,1)*[0:3])*1873);
    xStarSum = accumarray(xSumInds,xStar);
    profitFracMat = zeros(numel(lambdaX),ttMax);
    lambdaMat = zeros(numel(cityVectIndTrunc),ttMax);
    profitStar = xStar - (sum(uStar,1)' + (cV.*xStar));
    lambdaVect = zeros(numel(cityVectIndTrunc),ttMax);
    n = ones(size(T));
    U = uStar;

    jDMag = ones(size(jDStar));

    supplyFracTau = storagePhi;
    clear communityMatrixlongLouv xSumIndsAdd zFind tmFull zColIndFull timeMat

    %%%Run time loop
    %while ((xMag./xStarMag)<0.999)||(sum(shock)<(numCities.*numSectors))

    vI = T(:,:,1).*n(:,:,1);%mean cost per unit of input
    urbanToXInds = ones(numSectors,1)*[1:numCities]; urbanToXInds = urbanToXInds(:);

    productionAlpha = ones(size(xStar));
    inputAlpha = ones(size(xStar));

    P = ones(size(xStar));
    fdDistTotal = cStar;

    constructionDemands = zeros(size(lambdaX));
    constructionOutput = zeros(size(lambdaX));
    jDFull = jDFullStar;
    profitAlpha = (1./weeklyTau);
    runningProfitInds = repmat(([0:numCities-1]*numSectors),[numCities*numSectors 1]) + repmat([1:numSectors]',[numCities numCities]);
    cLAvg = ones(size(aStar));
    cLAvgVar = ones(size(aStar));
    totalIncomingDemandStar = accumarray(zStarRowInds,jDStar,[numCities*numSectors,1]);
    distanceMat = zeros(numSectors.*numCities);
    for i=1:numCities
        for j=1:numCities
            distanceMat([1:numSectors]+((i-1)*numSectors),[1:numSectors]+((j-1)*numSectors)) = time(i,j);
        end
    end

    clear distanceMat time
    constructionAllocationInds = repmat(([0:numCities-1]*numSectors)+constructionInd,[numCities*numSectors 1]);
    constructionAllocationInds = constructionAllocationInds(:);

    for tt = 1:ttMax
        %create impacts from the storm
        if tt==2
            currentDisasterInds = ((repmat(cityVectIndTrunc(find(urbanImpactMat(:,trackInd))),[1 numSectors])-1)*numSectors) + (ones(numel(find(urbanImpactMat(:,trackInd))),1)*[1:numSectors]);
            currentDisasterInds = currentDisasterInds';
            currentDisasterInds = currentDisasterInds(:);
            lambdaMag = repmat(wsMat(find(wsMat(:,trackInd)),trackInd),[1 numSectors])';
            lambdaMag = lambdaMag(:);
        else
            currentDisasterInds = 0;
            lambdaMag = 0;
        end

        fdiz = find(currentDisasterInds==0);
        currentDisasterInds(fdiz) = [];
        lambdaX(currentDisasterInds(:)) = lambdaX(currentDisasterInds(:)) + ((1-lambdaX(currentDisasterInds(:))).*lambdaMag(:));
        constructionDemands = lambdaX.*capital;

        %record lambda
        lambdaVect = accumarray(xSumInds(:),lambdaX,[numCities 1],@mean);
        lambdaMat(:,tt) = lambdaVect(cityVectIndTrunc);

        %incoming shipments from suppliers arrive
        I = T(:,:,1);
        vI = T(:,:,1).*n(:,:,1);%mean cost per unit of input
        vTotal = vS + sum(T.*n,3);
        sTotal = S + sum(T,3);
        unitaryCostTotal = (vTotal)./sTotal;
        unitaryCostTotal(sTotal==0) = 0;
        unitaryCostTotal(unitaryCostTotal<0) = 0;
        unitaryCost = vS./S;
        unitaryCost(S==0) = 0;
        unitaryCost(unitaryCost<0) = 0;
        vS = vI + vS;
        S = S + I;
        S(S<0) = 0;
        T(:,:,1:tStart-1) = T(:,:,2:end); T(:,:,end) = 0;
        n(:,:,1:tStart-1) = n(:,:,2:end); n(:,:,end) = 0;


        %calculate total demand for inputs for the next time
        %step from each purchaser. adjust based on the average
        %price of incoming goods
        deltaS = (sStar-S) + sum(tStar-T,3);
        totalConstructionDemands = zeros(size(lambdaX));
        fdDistTemp = cStar.*(unitaryCost.^-elasticity);
        fdDistTemp(unitaryCost==0) = 0;
        sRun = (sRun.*(1-(1./storagePhi))) + ((U + deltaS)./storagePhi);
        eU = U + deltaS + fdDistTemp;
        eU(eU<0) = 0;

        %calculate changes in the magnitude of demand
        dStar = cStar+uStar;
        jDMagTemp = eU./(dStar);
        jDMagTemp((dStar)==0) = 0;
        odaDenom = accumarray(xSumInds(:),sum(dStar,1)');
        odaDenomCorr = accumarray(xSumInds(:),sum(dStar,1)'.*(1-lambdaX));
        oda = accumarray(xSumInds(:),(sum(eU,1)'))./odaDenom;
        oda(odaDenom==0) = 1;
        odaCorr = accumarray(xSumInds(:),(sum(eU,1)'))./odaDenomCorr;
        odaCorr(odaDenomCorr==0) = 1;
        outgoingDemandAnomaly(:,tt) = oda(cityVectIndTrunc);
        outgoingDemandAnomalyCorr(:,tt) = odaCorr(cityVectIndTrunc);
        jDMag = jDMagTemp(uStarInds);

        %calculate demand that can be satisfied by suppliers in
        %next time step w/ potential increase in alpha
        inputFrac = S./uStar;
        uColInds = repmat([1:numCities*numSectors],[numSectors,1]);
        inff = inputFrac; inff(uStar==0) = 0;
        inputFrac(uStar==0) = max(inff(:,uColInds(uStar==0)),[],1)';
        inputCapacity = min(inputFrac,[],1)';
        prospectiveDemandFill = min([(productionAlpha.*(1-lambdaX)) inputAlpha.*inputCapacity],[],2);
        constructionDemandFrac = totalConstructionDemands./xStar;
        constructionDemandFrac(xStar==0) = 0;
        prospectiveDemandFill = prospectiveDemandFill - constructionDemandFrac;
        prospectiveDemandFill(prospectiveDemandFill<0) = 0;
        prospectiveDemandFillZ = prospectiveDemandFill(zStarRowInds);
        pElasticity = 1 - (demandElasticity.*(P-1));
        pElasticity(pElasticity<0) = 0;
        pElasticityZ = pElasticity(zStarRowInds);
        %distribute demand to suppliers
        preferredjDGamma = ones(size(jDStar));
        supplierPricePreferred = zeros(size(jDStar));
        PZ = P(zStarRowInds);
        for i=1:numSectors*numCities
            for j=1:numSectors
                supplierInds = jStruct(i).sector(j).inds;
                currentjD = jDStar(supplierInds);
                suppliers = currentjD.*prospectiveDemandFillZ(supplierInds).*pElasticityZ(supplierInds);
                if sum(suppliers)>0
                    preferredConfiguration = suppliers./sum(suppliers);
                    jDTemp = preferredConfiguration.*dStar(j,i);
                    preferredjDGamma(supplierInds) = jDTemp./currentjD;
                    supplierPricePreferred(supplierInds) = 1-(demandElasticity.*(sum(PZ(supplierInds).*preferredConfiguration)-1));
                end
            end
        end
        %calculate an overall jD
        jDGamma = jDGamma + ((preferredjDGamma-jDGamma).*(1-exp(-1./weeklyTau)));
        jD = jDGamma.*jDMag.*jDStar.*supplierPricePreferred;

        %record changes in jD for each urban area associated with
        %those urban areas directly affected
        purchaserIncomingDemand = accumarray(zStarRowInds,jD);

        %calculate input costs
        cL = sum(unitaryCostTotal.*aStar,1)';%commodity cost of inputs
        cL(xStarCol==0) = 0;

        %adjust alpha to accomodate changes in demand that are
        %profitable
        totalIncomingDemand = purchaserIncomingDemand;
        targetAlpha = (totalIncomingDemand)./xStar;
        targetAlpha(xStar==0) = 1;
        targetAlpha(targetAlpha<1) = 1; targetAlpha(targetAlpha>alphaMax) = alphaMax;
        idaDenom = accumarray(xSumInds,xStar);
        idaDenomCorr = accumarray(xSumInds,xStar.*(1-lambdaX));
        ida = accumarray(xSumInds,totalIncomingDemand)./idaDenom;
        ida(idaDenom==0) = 1;
        idaCorr = accumarray(xSumInds,totalIncomingDemand)./idaDenomCorr;
        idaCorr(idaDenomCorr==0) = 1;
        incomingDemandAnomaly(:,tt) = ida(cityVectIndTrunc);
        incomingDemandAnomalyCorr(:,tt) = idaCorr(cityVectIndTrunc);
        profitCalc = P - ((cL./inputAlpha)+cV);
        targetAlphaProfit = zeros(size(targetAlpha));
        for i=1:numel(xStar)
            if profitCalc(i)<0
                targetAlphaProfit(i) = 1;
            else
                profitTemp = (P(i) - (cL(i)./inputAlpha(i))) - (cVSub(i) + (cVCost(i).*linspace(1,targetAlpha(i),100).*inputAlpha(i)));
                targetAlphaProfit(i) = 1 + ((sum(profitTemp>0)./100)*(targetAlpha(i)-1));
            end
        end

        productionAlpha = productionAlpha + ((targetAlphaProfit - productionAlpha).*(1-exp(-1./overallTau)));
        actualProfitRate = round(P - ((cL./inputAlpha) + (cV.*productionAlpha)),3);
        %inputTarget = targetAlphaProfit./inputCapacity;
        %inputTarget(inputCapacity==0) = 1; inputTarget(inputTarget<1) = 1; inputTarget(inputTarget>alphaMax) = alphaMax;
        profitFrac = profitFrac + (((actualProfitRate>=0)-profitFrac).*(1-exp(-1./overallTau)));

        %produce output based on demand
        [cc limitationVect1] = min([(1-lambdaX).*productionAlpha inputAlpha.*inputCapacity],[],2);
        currentCapacity = cc.*xStar.*profitFrac;
        zOut = zeros(size(jDStar));
        X = zeros(size(xStar));
        revenue = zeros(size(xStar));
        producedFrac = currentCapacity./totalIncomingDemand;
        producedFrac(totalIncomingDemand==0) = 1;
        demandFrac = totalIncomingDemand./totalIncomingDemandStar;
        limitationArray = ([(1-lambdaX).*productionAlpha inputAlpha.*inputCapacity profitFrac demandFrac]).*repmat(xStar,1,4);
        [limitationVect] = (reshape(accumarray(limInds(:),limitationArray(:)),1873,1,4));

        constructionOutput = zeros(size(constructionDemands));
        constructionDemandFrac = zeros(size(totalIncomingDemand));
        constructionDemandFrac(totalIncomingDemand==0) = 0;
        for i=1:numSectors*numCities
            if xStar(i)>0
                zRowInds = zRowStruct(i).inds;
                if producedFrac(i)<1
                    X(i) = totalIncomingDemand(i).*producedFrac(i);
                    revenue(i) = X(i).*P(i);
                    constructionOutput(i) = X(i).*constructionDemandFrac(i);
                    zOut(zRowInds) = (jD(zRowInds)./sum(jD(zRowInds))).*X(i).*(1-constructionDemandFrac(i));
                else
                    constructionOutput(i) = totalConstructionDemands(i);
                    zOut(zRowInds) = jD(zRowInds);
                    X(i) = totalIncomingDemand(i);
                    revenue(i) = X(i).*P(i);
                end
            end
        end

        xDepressed = xStar-X;
        xDepressed(xDepressed<0) = 0;
        tid = totalIncomingDemand;
        tid(tid>xStar) = xStar(tid>xStar);

        %load goods and value onto transit matrix
        zTrans = reshape(accumarray(tInds,zOut,[numSectors*numSectors*numCities*maxDistance,1]),numSectors,numSectors*numCities,maxDistance);
        nTransNum = reshape(accumarray(tInds,P(zStarRowInds).*zOut,[numSectors*numSectors*numCities*maxDistance,1]),numSectors,numSectors*numCities,maxDistance);
        nTrans = nTransNum./zTrans;
        nTrans(zTrans==0) = 0;
        n = ((nTrans.*zTrans) + (n.*T))./(zTrans+T);
        n((zTrans+T)==0) = 1;
        T = T + zTrans;

        %adjust storage for inputs to production and to satisfy
        %FD
        U = uStar.*repmat((X./(xStar.*inputAlpha))',[numSectors,1]);
        U(:,xStar==0) = 0;
        unitaryCost = (vS)./S;
        unitaryCost(S==0) = 0;
        unitaryCost(unitaryCost<0) = 0;
        vS = vS - (U.*unitaryCost);
        S = S - U;
        S(S<0) = 0;%account for rounding errors
        vS(S==0) = 0;

        %distribute to FD
        fdDistTotal = cStar.*((unitaryCost./1).^-elasticity);
        fdDistTotal(unitaryCost==0) = 0;

        fdDistSum = fdDistTotal;
        fractionFDServed = (S./storagePhi)./(fdDistSum);
        fractionFDServed((fdDistSum)==0) = 0;
        fractionFDServed(fractionFDServed>1) =1;
        vS = vS - (fdDistSum.*unitaryCost.*fractionFDServed);
        S = S - (fdDistSum.*fractionFDServed);
        fdTotalFrac = fdDistTotal./(fdDistSum);
        fdTotalFrac(fdDistSum==0) = 0;
        fdDistTotal = fdDistTotal.*fdTotalFrac;

        %adjust construction matrix
        constructionOutputToFD = constructionOutput;
        constructionDemandFrac = (1/rebuildTau);
        constructionDemandFrac(totalConstructionDemands==0) = 0;
        constructionDemandFrac(constructionDemandFrac<0) = 0;
        constructionDemands = constructionDemands.*(1-(1/overallTau));
        lambdaX = constructionDemands./capital;
        lambdaX(capital==0) = 0;

        %adjust price for the next time step
        P = (1 + (gammaP.*(totalIncomingDemand - X)./X));
        P(X==0) = 1+gammaP;

        %save some state variables for plotting
        xsix = accumarray(xSumInds(:),X);
        totalXVect(:,tt) = xsix(cityVectIndTrunc);
        limitationMat(:,tt,:) = limitationVect(cityVectIndTrunc,:,:);
        tt
    end

    save(sprintf('%srunOutputStorm%dParameter%dTest%d.mat',outputPath,trackInd,paramNum,paramVal),'lambdaMat','totalXVect','lambdaVect','xStar');

end
