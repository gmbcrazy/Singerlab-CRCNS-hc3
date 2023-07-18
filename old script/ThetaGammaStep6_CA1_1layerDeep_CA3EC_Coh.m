        

clear all
%%%%%%%%load theta filter parameters created by ThetaFilterAll_BKdata.m 
%%%%%%%%load theta filter parameters created by ThetaFilterAll_BKdata.m 
rawFolder='Y:\rozell-singer-exchange\CRCNSdata\hc-3\dataAll\';
proccessFolder='Y:\rozell-singer-exchange\CRCNSdata\hc-3\ProccessedDataAll\';
close all
ProcessGroup=dir([proccessFolder 'CA1*']);
samprate=1250;%%sampling rate
WaveParam.ntw=11;      %smoothing time-window 11
WaveParam.nsw=3;        %smoothing fre-window 3
WaveParam.Range=0.2;      %AlignedWindow Size 0.2s
WaveParam.smoothW=10;   %smoothing Power
WaveParam.DownSample=2; %DownSampling
WaveParam.wname='morl'; %Wavelet Name
WaveParam.Samplingrate=samprate;
WaveParam.Freq=[20:2:180]; %Frequeny band interested
WaveParam.Zscore=1;

PhaseNum=20;
PhaseStep=2*pi/PhaseNum;
Param.PhaseBin=-pi:PhaseStep:pi;
Cnum=4;

 %%%%%%define ploting colors for different cluster
 ComColor= [10,103,155;176,34,242;78,211,34;191,126,0;178,51,51;51,204,178]/255;
 %%%%%%define ploting colors for different cluster
%%%%%%%%%%parameter for visulization of aveaging sample
P.Ytick=[20:40:180];
P.Yticklabel=[20:40:180];
P.Xtick=[-pi 0 pi];
P.Xticklabel={'-\pi','0','\pi'};
P.Ylabel='Frequency Hz';
P.Xlabel='Theta Phase rad';
P.CbarYLabel='Coherence';
P.Clim=[0 0.6];

load('Y:\rozell-singer-exchange\CRCNSdata\hc-3\AnalysisAll\MatchResult\AnimalBehaviorAll.mat');
for i=1:length(ProcessGroup)
    temp1=[ProcessGroup(i).folder '\' ProcessGroup(i).name '\'];
    FileL1=dir(temp1);
    FileL1(1:2)=[];
    Invalid=[];
    for itop=1:length(FileL1)
        if ~isempty(strfind(FileL1(itop).name,'Used'))
           Invalid=[Invalid;itop];
        elseif FileL1(itop).isdir~=1
           Invalid=[Invalid;itop]; 
        else    
        end
    end
    FileL1(Invalid)=[];
    clear itop
    
    

    for itop=1:length(FileL1) 
        FileL2=dir([FileL1(itop).folder '\' FileL1(itop).name]);
        FileL2(1:2)=[];
        N1=[];
        for itemp=1:length(FileL2)
            if FileL2(itemp).isdir
               N1=[N1;itemp]; 
            end
        end
        FileL2=FileL2(N1);

        load([FileL2(1).folder '\PyrLayerInfo1.mat']);
        load([FileL2(1).folder '\BestCA3ECChan_UniqueChPerShank.mat']);
        
        CA1Chan=AnimalGroup(i).Animal(itop).PyrChan;
        CA1Shank=PyrShank(:,2);
        numCA1=length(CA1Chan);
        
        


        for iS=1:length(FileL2)
        PathTemp=[FileL2(iS).folder '\' FileL2(iS).name '\'];
% %         load([PathTemp 'thetaEpoch.mat']);
        load([PathTemp 'SleepState.mat']);
        if ~isempty(strfind(ProcessGroup(i).name,'sleep'))
           Period=SleepState.ints.REMstate';
        elseif ~isempty(strfind(ProcessGroup(i).name,'linear'))
           Period=SleepState.ints.WAKEstate';
        else
           error('NoPeriod Detected');
        end

            
            
            for iCA1=1:numCA1
                temp=[FileL2(1).folder '\' AnimalGroup(i).Animal(itop).Behavior(iCA1).Period(iS).SessionName '\'];
                kCA1=['theta_Ch' num2str(CA1Chan(iCA1)) 'Sh' num2str(CA1Shank(iCA1)) 'CA1.mat'];
                load([temp kCA1]);
                CA1theta=lfpTheta; 
                clear lfpTheta
                CA1phase=phase;
                clear phase
                
                %%%%%One Layer Deeper than Pyrlayer for Coh-Calculation
                kCA1=['Ch' num2str(CA1Chan(iCA1)+1) 'Sh' num2str(CA1Shank(iCA1)) 'CA1.mat'];
                %%%%%
                load([temp kCA1]);
                 CA1lfp=eeg;
                 clear eeg

                for iCAX=1:length(ShankCh)
                    if ShRegion(iCAX)>0
                        key1=['CA' num2str(ShRegion(iCAX)) '.mat'];
                    elseif ShRegion(iCAX)<0
                        key1=['EC' num2str(abs(ShRegion(iCAX))) '.mat'];
                    else      
                    end
                    
                    load([temp 'Ch' num2str(ShankCh(iCAX)) 'Sh' num2str(ShankSave(iCAX)) key1]);
                    CAXlfp=eeg;
                    clear eeg;
                    Param.SavePath=PathTemp;
                    Param.SaveShow=[ProcessGroup(i).name '_Coh1LayerDeepCh' num2str(ShankCh(iCAX)) 'PhaseCh' num2str(CA1Chan(iCA1))];

                    [~,~] = LU_CohThetaExtract(CA1lfp,CAXlfp,CA1theta,CA1phase,samprate,Period,WaveParam,Param);
                    KeyWordSave=[kCA1 'and' key1];
                    KeyWordAll(i).Animal(itop).KeyWord{iS,iCA1,iCAX}=KeyWordSave;
                end
            end
            
        end
    end

end


for i=1:length(ProcessGroup)
    temp1=[ProcessGroup(i).folder '\' ProcessGroup(i).name '\'];
    FileL1=dir(temp1);
    FileL1(1:2)=[];
    Invalid=[];
    for itop=1:length(FileL1)
        if ~isempty(strfind(FileL1(itop).name,'Used'))
           Invalid=[Invalid;itop];
        elseif FileL1(itop).isdir~=1
           Invalid=[Invalid;itop]; 
        else    
        end
    end
    FileL1(Invalid)=[];
    clear itop
    
    

    for itop=1:length(FileL1) 
        FileL2=dir([FileL1(itop).folder '\' FileL1(itop).name]);
        FileL2(1:2)=[];
        N1=[];
        for itemp=1:length(FileL2)
            if FileL2(itemp).isdir
               N1=[N1;itemp]; 
            end
        end
        FileL2=FileL2(N1);

        load([FileL2(1).folder '\PyrLayerInfo1.mat']);
        load([FileL2(1).folder '\BestCA3ECChan_UniqueChPerShank.mat']);
        
        CA1Chan=AnimalGroup(i).Animal(itop).PyrChan;
        CA1Shank=PyrShank(:,2);
        numCA1=length(CA1Chan);
        
        CAXChan=ShankCh;
        numCAX=length(CAXChan);
        clear CAXName
        for iCAX=1:numCAX
            if ShRegion(iCAX)<0
               CAXName{iCAX}=['EC' num2str(abs(ShRegion(iCAX)))];
            else
               CAXName{iCAX}=['CA3' num2str(abs(ShRegion(iCAX)))];
            end
            
        end
        
        for iCA1=1:numCA1
           for iCAX=1:numCAX
               sampleA=[];
               for iS=1:length(FileL2)
                   PathTemp=[FileL2(iS).folder '\' FileL2(iS).name '\'];
                   TempFile=[ProcessGroup(i).name '_Coh1LayerDeepCh' num2str(CAXChan(iCAX)) ...
                       'PhaseCh' num2str(CA1Chan(iCA1))];
                   load([PathTemp TempFile]);
                   sampleA=[sampleA sample];
               end
               sampleN=size(sampleA,2);
               TempSS=AnimalGroup(i).Animal(itop).Behavior(iCA1).SortedTsState;
               TempCom=AnimalGroup(i).Animal(itop).Behavior(iCA1).SortedCom;
               PhasePlot=AnimalGroup(i).Animal(itop).PhasePlot;
               DataCom=[];
               if length(TempSS)~=sampleN
                   error('NotMatched')
               else
                   for iCom=1:Cnum
                       DataCom(:,:,iCom)=reshape(mean(sampleA(:,TempSS==iCom),2),Div(1),Div(2));
                   end
        figure;
        MultiMatrix2DPlot(abs(DataCom),1:Cnum,PhasePlot,Fplot,P);
        
%         imagesc(PhasePlot,Fplot,DataCom(:,:,3)');axis xy
        
        colormap(jet);
        papersizePX=[0 0 20 7];
        set(gcf, 'PaperUnits', 'centimeters');
        set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
        SaveName=[FileL2(1).folder '\' CAXName{iCAX} 'Coh1LayerDeepCh' num2str(CAXChan(iCAX)) 'PhaseCh' num2str(CA1Chan(iCA1))];
        saveas(gcf,[SaveName '.eps'],'epsc'); 
        saveas(gcf,[SaveName],'tiff'); 
        close all
                  AnimalGroup(i).Animal(itop).Behavior(iCA1).CAX(iCAX).Com=DataCom;
               end
           end
        end
    end

end

save('Y:\rozell-singer-exchange\CRCNSdata\hc-3\AnalysisAll\MatchResult\AnimalCA1Coh1LayerDeep_CA3EC.mat','KeyWordAll','AnimalGroup');

load('Y:\rozell-singer-exchange\CRCNSdata\hc-3\AnalysisAll\MatchResult\AnimalCA1Coh_CA3EC.mat','KeyWordAll','AnimalGroup');

for i=1:2
    rAll=[];

    for itop=1:length(AnimalGroup(i).Animal) 

%         load([FileL2(1).folder '\PyrLayerInfo1.mat']);
%         load([FileL2(1).folder '\BestCA3ECChan_UniqueChPerShank.mat']);
%         
%         CA1Chan=AnimalGroup(i).Animal(itop).PyrChan;
%         CA1Shank=PyrShank(:,2);
        numCA1=length(AnimalGroup(i).Animal(itop).Behavior);
        
%         CAXChan=ShankCh;
        numCAX=length(AnimalGroup(i).Animal(itop).Behavior(1).CAX);
        
        for iCA1=1:numCA1
            Dim=size(AnimalGroup(i).Animal(itop).Behavior(iCA1).SortedCom);
            tempCA1=reshape(AnimalGroup(i).Animal(itop).Behavior(iCA1).SortedCom(:,1:end,:),Dim(1)*(Dim(2)),Dim(3));
            Cov=reshape(nanmean(AnimalGroup(i).Animal(itop).Behavior(iCA1).SortedCom(:,1:end,:),3),Dim(1)*(Dim(2)),1);

            rCA=zeros(Cnum,Cnum);
%                 figure;
%                 MultiMatrix2DPlot(AnimalGroup(i).Animal(itop).Behavior(iCA1).SortedCom,1:Cnum,PhasePlot,Fplot);
%                 colormap(hot)
            CountCAX=0;
            for iCAX=1:numCAX     
                temp=AnimalGroup(i).Animal(itop).Behavior(iCA1).CAX(iCAX).Com;
                Dim=size(temp);
                tempCAX=reshape(AnimalGroup(i).Animal(itop).Behavior(iCA1).CAX(iCAX).Com(:,1:end,:),Dim(1)*(Dim(2)),Dim(3));
                r=partialcorr(tempCA1,abs(tempCAX),Cov);
                rCA=rCA+r;
                CountCAX=CountCAX+1;
%                 figure
%                 MultiMatrix2DPlot(abs(AnimalGroup(i).Animal(itop).Behavior(iCA1).CAX(iCAX).Com),1:Cnum,PhasePlot,Fplot);
%                 colormap(jet)

            end
            
            if CountCAX>0
            rCA=rCA/CountCAX;
            rAll=cat(3,rAll,rCA);
            end

        end

    end
    rCohPower_CA3CA1(i).Animal=rAll;

% %     papersizePX=[0 0 8*numCA1 8];
% %     set(gcf, 'PaperUnits', 'centimeters');
% %     set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
% %     SaveName=[FileL2(1).folder '\' 'CA1TriggeredResult'];
% % % %     saveas(gcf,[SaveName '.eps'],'epsc'); 
% %     saveas(gcf,[SaveName],'tiff'); 

end
for i=3:4
    rAll=[];

    for itop=1:length(AnimalGroup(i).Animal) 

%         load([FileL2(1).folder '\PyrLayerInfo1.mat']);
%         load([FileL2(1).folder '\BestCA3ECChan_UniqueChPerShank.mat']);
        
%         CA1Chan=AnimalGroup(i).Animal(itop).PyrChan;
%         CA1Shank=PyrShank(:,2);
        numCA1=length(AnimalGroup(i).Animal(itop).Behavior); 
%         CAXChan=ShankCh;
        numCAX=length(AnimalGroup(i).Animal(itop).Behavior(1).CAX);
        
        for iCA1=1:numCA1
            Dim=size(AnimalGroup(i).Animal(itop).Behavior(iCA1).SortedCom);
            tempCA1=reshape(AnimalGroup(i).Animal(itop).Behavior(iCA1).SortedCom(:,1:end,:),Dim(1)*(Dim(2)),Dim(3));
            rCA=zeros(Cnum,Cnum);
%                 figure;
%                 MultiMatrix2DPlot(AnimalGroup(i).Animal(itop).Behavior(iCA1).SortedCom,1:Cnum,PhasePlot,Fplot);
%                 colormap(hot)
            CountCAX=0;
            for iCAX=1:numCAX    
                if isempty(strfind(KeyWordAll(i).Animal(itop).KeyWord{1,iCA1,iCAX},'EC2'))  %%%%%%%%EC2 
                   continue;
                end
                CountCAX=CountCAX+1;
                
                temp=AnimalGroup(i).Animal(itop).Behavior(iCA1).CAX(iCAX).Com;
                Dim=size(temp);
                tempCAX=reshape(AnimalGroup(i).Animal(itop).Behavior(iCA1).CAX(iCAX).Com(:,1:end,:),Dim(1)*(Dim(2)),Dim(3));
                r=corr(tempCA1,abs(tempCAX));
                rCA=rCA+r;
%                 figure
%                 MultiMatrix2DPlot(abs(AnimalGroup(i).Animal(itop).Behavior(iCA1).CAX(iCAX).Com),1:Cnum,PhasePlot,Fplot);
%                 colormap(jet)

            end
            if CountCAX>0
            rCA=rCA/CountCAX;
            rAll=cat(3,rAll,rCA);
            end
        end

    end
    rCohPower_EC2CA1(i-2).Animal=rAll;

% %     papersizePX=[0 0 8*numCA1 8];
% %     set(gcf, 'PaperUnits', 'centimeters');
% %     set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
% %     SaveName=[FileL2(1).folder '\' 'CA1TriggeredResult'];
% % % %     saveas(gcf,[SaveName '.eps'],'epsc'); 
% %     saveas(gcf,[SaveName],'tiff'); 

end
for i=3:4
    rAll=[];

    for itop=1:length(AnimalGroup(i).Animal) 

%         load([FileL2(1).folder '\PyrLayerInfo1.mat']);
%         load([FileL2(1).folder '\BestCA3ECChan_UniqueChPerShank.mat']);
%         
%         CA1Chan=AnimalGroup(i).Animal(itop).PyrChan;
%         CA1Shank=PyrShank(:,2);
        numCA1=length(AnimalGroup(i).Animal(itop).Behavior);
        
%         CAXChan=ShankCh;
        numCAX=length(AnimalGroup(i).Animal(itop).Behavior(1).CAX);
        
        for iCA1=1:numCA1
            Dim=size(AnimalGroup(i).Animal(itop).Behavior(iCA1).SortedCom);
            tempCA1=reshape(AnimalGroup(i).Animal(itop).Behavior(iCA1).SortedCom(:,1:end,:),Dim(1)*(Dim(2)),Dim(3));
            rCA=zeros(Cnum,Cnum);
%                 figure;
%                 MultiMatrix2DPlot(AnimalGroup(i).Animal(itop).Behavior(iCA1).SortedCom,1:Cnum,PhasePlot,Fplot);
%                 colormap(hot)
            CountCAX=0;
            for iCAX=1:numCAX    
                if isempty(strfind(KeyWordAll(i).Animal(itop).KeyWord{1,iCA1,iCAX},'EC3'))  %%%%%%%%EC2 
                   continue;
                end
                CountCAX=CountCAX+1;
                
                temp=AnimalGroup(i).Animal(itop).Behavior(iCA1).CAX(iCAX).Com;
                Dim=size(temp);
                tempCAX=reshape(AnimalGroup(i).Animal(itop).Behavior(iCA1).CAX(iCAX).Com(:,1:end,:),Dim(1)*(Dim(2)),Dim(3));
                r=corr(tempCA1,abs(tempCAX));
                rCA=rCA+r;
%                 figure
%                 MultiMatrix2DPlot(abs(AnimalGroup(i).Animal(itop).Behavior(iCA1).CAX(iCAX).Com),1:Cnum,PhasePlot,Fplot);
%                 colormap(jet)

            end
            if CountCAX>0
            rCA=rCA/CountCAX;
            rAll=cat(3,rAll,rCA);
            end

        end

    end
    rCohPower_EC3CA1(i-2).Animal=rAll;

% %     papersizePX=[0 0 8*numCA1 8];
% %     set(gcf, 'PaperUnits', 'centimeters');
% %     set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
% %     SaveName=[FileL2(1).folder '\' 'CA1TriggeredResult'];
% % % %     saveas(gcf,[SaveName '.eps'],'epsc'); 
% %     saveas(gcf,[SaveName],'tiff'); 

end
for i=3:4
    rAll=[];

    for itop=1:length(AnimalGroup(i).Animal) 

% %         load([FileL2(1).folder '\PyrLayerInfo1.mat']);
% %         load([FileL2(1).folder '\BestCA3ECChan_UniqueChPerShank.mat']);
        
%         CA1Chan=AnimalGroup(i).Animal(itop).PyrChan;
%         CA1Shank=PyrShank(:,2);
        numCA1=length(AnimalGroup(i).Animal(itop).Behavior);
%         CAXChan=ShankCh;
        numCAX=length(AnimalGroup(i).Animal(itop).Behavior(1).CAX);
        
        for iCA1=1:numCA1
            Dim=size(AnimalGroup(i).Animal(itop).Behavior(iCA1).SortedCom);
            tempCA1=reshape(AnimalGroup(i).Animal(itop).Behavior(iCA1).SortedCom(:,1:end,:),Dim(1)*(Dim(2)),Dim(3));
            rCA=zeros(Cnum,Cnum);
%                 figure;
%                 MultiMatrix2DPlot(AnimalGroup(i).Animal(itop).Behavior(iCA1).SortedCom,1:Cnum,PhasePlot,Fplot);
%                 colormap(hot)
            CountCAX=0;
            for iCAX=1:numCAX    
                if isempty(strfind(KeyWordAll(i).Animal(itop).KeyWord{1,iCA1,iCAX},'EC4'))  %%%%%%%%EC2 
                   continue;
                end
                CountCAX=CountCAX+1;
                
                temp=AnimalGroup(i).Animal(itop).Behavior(iCA1).CAX(iCAX).Com;
                Dim=size(temp);
                tempCAX=reshape(AnimalGroup(i).Animal(itop).Behavior(iCA1).CAX(iCAX).Com(:,1:end,:),Dim(1)*(Dim(2)),Dim(3));
                r=corr(tempCA1,abs(tempCAX));
                rCA=rCA+r;
%                 figure
%                 MultiMatrix2DPlot(abs(AnimalGroup(i).Animal(itop).Behavior(iCA1).CAX(iCAX).Com),1:Cnum,PhasePlot,Fplot);
%                 colormap(jet)

            end
            if CountCAX>0
            rCA=rCA/CountCAX;
            rAll=cat(3,rAll,rCA);
            end

        end

    end
    rCohPower_EC4CA1(i-2).Animal=rAll;

% %     papersizePX=[0 0 8*numCA1 8];
% %     set(gcf, 'PaperUnits', 'centimeters');
% %     set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
% %     SaveName=[FileL2(1).folder '\' 'CA1TriggeredResult'];
% % % %     saveas(gcf,[SaveName '.eps'],'epsc'); 
% %     saveas(gcf,[SaveName],'tiff'); 

end
for i=3:4
    rAll=[];

    for itop=1:length(AnimalGroup(i).Animal) 

% %         load([FileL2(1).folder '\PyrLayerInfo1.mat']);
% %         load([FileL2(1).folder '\BestCA3ECChan_UniqueChPerShank.mat']);
        
%         CA1Chan=AnimalGroup(i).Animal(itop).PyrChan;
%         CA1Shank=PyrShank(:,2);
        numCA1=length(AnimalGroup(i).Animal(itop).Behavior);
%         CAXChan=ShankCh;
        numCAX=length(AnimalGroup(i).Animal(itop).Behavior(1).CAX);
        
        for iCA1=1:numCA1
            Dim=size(AnimalGroup(i).Animal(itop).Behavior(iCA1).SortedCom);
            tempCA1=reshape(AnimalGroup(i).Animal(itop).Behavior(iCA1).SortedCom(:,1:end,:),Dim(1)*(Dim(2)),Dim(3));
            rCA=zeros(Cnum,Cnum);
%                 figure;
%                 MultiMatrix2DPlot(AnimalGroup(i).Animal(itop).Behavior(iCA1).SortedCom,1:Cnum,PhasePlot,Fplot);
%                 colormap(hot)
            CountCAX=0;
            for iCAX=1:numCAX    
                if isempty(strfind(KeyWordAll(i).Animal(itop).KeyWord{1,iCA1,iCAX},'EC5'))  %%%%%%%%EC2 
                   continue;
                end
                CountCAX=CountCAX+1;
                
                temp=AnimalGroup(i).Animal(itop).Behavior(iCA1).CAX(iCAX).Com;
                Dim=size(temp);
                tempCAX=reshape(AnimalGroup(i).Animal(itop).Behavior(iCA1).CAX(iCAX).Com(:,1:end,:),Dim(1)*(Dim(2)),Dim(3));
                r=corr(tempCA1,abs(tempCAX));
                rCA=rCA+r;
%                 figure
%                 MultiMatrix2DPlot(abs(AnimalGroup(i).Animal(itop).Behavior(iCA1).CAX(iCAX).Com),1:Cnum,PhasePlot,Fplot);
%                 colormap(jet)

            end
            if CountCAX>0
            rCA=rCA/CountCAX;
            rAll=cat(3,rAll,rCA);
            end

        end

    end
    rCohPower_EC5CA1(i-2).Animal=rAll;

% %     papersizePX=[0 0 8*numCA1 8];
% %     set(gcf, 'PaperUnits', 'centimeters');
% %     set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
% %     SaveName=[FileL2(1).folder '\' 'CA1TriggeredResult'];
% % % %     saveas(gcf,[SaveName '.eps'],'epsc'); 
% %     saveas(gcf,[SaveName],'tiff'); 

end

% temp1=rCohPower_CA3CA1(1);
Data(1,1).Comparison={};
Data(1,1).Name='LinearCA3CA1';
Data(1,2).Comparison={};
Data(1,2).Name='LinearEC2CA1';
Data(1,3).Comparison={};
Data(1,3).Name='LinearEC3CA1';
Data(1,4).Comparison={};
Data(1,4).Name='LinearEC4CA1';
Data(1,5).Comparison={};
Data(1,5).Name='LinearEC5CA1';
% temp1=rCohPower_CA3CA1(1);
Data(2,1).Comparison={};
Data(2,1).Name='SleepCA3CA1';
Data(2,2).Comparison={};
Data(2,2).Name='SleepEC2CA1';
Data(2,3).Comparison={};
Data(2,3).Name='SleepEC3CA1';
Data(2,4).Comparison={};
Data(2,4).Name='SleepEC4CA1';
Data(2,5).Comparison={};
Data(2,5).Name='SleepEC5CA1';


for i=1:size(Data,1)
for iPower=1:Cnum
    for iCoh=1:Cnum
%         temp2{end+1}=rCohPower_CA3CA1(1).Animal(iPower,iCoh,:);
        Data(i,1).Comparison{end+1}=rCohPower_CA3CA1(i).Animal(iPower,iCoh,:);
        Data(i,2).Comparison{end+1}=rCohPower_EC2CA1(i).Animal(iPower,iCoh,:);
        Data(i,3).Comparison{end+1}=rCohPower_EC3CA1(i).Animal(iPower,iCoh,:);
        Data(i,4).Comparison{end+1}=rCohPower_EC4CA1(i).Animal(iPower,iCoh,:);
        Data(i,5).Comparison{end+1}=rCohPower_EC5CA1(i).Animal(iPower,iCoh,:);
    end
end
end
Datatype=2;
x=[1:4 6:9 11:14 16:19];
% stats=ErrorBarPlotLU(x,temp2,1,repmat(ComColor(1:4,:),4,1),2)

Datatype=0;

   GroupPair.CorrName='fdr';
   GroupPair.Q=0.1;
   GroupPair.Pair=[];
   GroupPair.SignY=0.6;
GroupPair.Plot=0;
figure;
for i=1:size(Data,1)
    for j=1:size(Data,2)
        subplotLU(size(Data,1),size(Data,2),i,j)
        stats=ErrorBarPlotLU(x,Data(i,j).Comparison,1,repmat(ComColor(1:4,:),4,1),2,1,[],GroupPair);
    end    
end




