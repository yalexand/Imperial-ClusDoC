function [CellData_out, DensityROI] = DoCHandler_YA(ROICoordinates, CellData_in, Lr_rad, Rmax, Step, Chan1Color, Chan2Color, Outputfolder, NDatacolumns)

    % Handler function for data arrangement and plotting of DoC calculations
    % Formerly Main_Fun_DoC_GUIV2.m
    % Follows idea of:
    % Colocalization Malkusch (Histochem Cell Biol)
    % http://www.ncbi.nlm.nih.gov/pubmed/22086768
    % Moved here to avoid conflicts with private functions and data
    % handling.
    % Away from button callback for clarity.

    CellData_out = CellData_in;
    
    Data_DoC = {[]};    % cell(max(cellfun(@length, ROICoordinates)), length(CellData));
    DensityROI = {[]};  % cell(max(cellfun(@length, ROICoordinates)), length(CellData));

        for roiIter = 1:numel(ROICoordinates) % ROI number

                    roi = ROICoordinates{roiIter};                       
                    x_roi=roi(:,1);
                    y_roi=roi(:,2);
                        x=CellData_in(:,5);
                        y=CellData_in(:,6);
                    whichPointsInROI = x>=min(x_roi) & x<=max(x_roi) & y>=min(y_roi) & y<=max(y_roi);                            
                    dataCropped = CellData_in(whichPointsInROI,:);                   

            if ~isempty(dataCropped)
                
                dataCropped(isnan(dataCropped(:,12)),:)=[];

                % Send single ROI of data to DoC calculation function
                [ Data_DegColoc1, SizeROI1 ] = DoCCalc( dataCropped, Lr_rad, Rmax, Step, roi );

                CA1 = Data_DegColoc1.DoC((Data_DegColoc1.Ch == 1) & (Data_DegColoc1.Lr_rAboveThresh == 1)); % Ch1 -> Ch2
                CA2 = Data_DegColoc1.DoC((Data_DegColoc1.Ch == 2) & (Data_DegColoc1.Lr_rAboveThresh == 1)); % Ch2 -> Ch1

                handles.handles.DoCFigPerROI = figure('color', [1 1 1], 'inverthardcopy', 'off');
                                
                handles.handles.DoCAxPerROI(1) = subplot(2,1,1);
                histHand = histogram(handles.handles.DoCAxPerROI(1), CA1, 100);
                histHand.FaceColor = Chan1Color;
                histHand.EdgeColor = rgb(52, 73, 94);
                set(handles.handles.DoCAxPerROI(1), 'XLim', [-1 1]);
                xlabel(handles.handles.DoCAxPerROI(1), 'DoC Score Ch1 -> Ch2', 'Fontsize', 20);
                ylabel(handles.handles.DoCAxPerROI(1), 'Frequency','FontSize',20);
                set(handles.handles.DoCAxPerROI(1),'FontSize',20)

                handles.handles.DoCAxPerROI(2) = subplot(2,1,2);
                histHand = histogram(handles.handles.DoCAxPerROI(2), CA2, 100);
                histHand.FaceColor = Chan2Color;
                histHand.EdgeColor = rgb(52, 73, 94);
                set(handles.handles.DoCAxPerROI(2), 'XLim', [-1 1]);
                xlabel(handles.handles.DoCAxPerROI(2), 'DoC Score Ch2 -> Ch1', 'Fontsize', 20);
                ylabel(handles.handles.DoCAxPerROI(2), 'Frequency','FontSize',20);
                set(handles.handles.DoCAxPerROI(2),'FontSize',20)

                drawnow;

                % Save the figure
                Name = sprintf('Region_%d_Hist', roiIter);
                print(fullfile(Outputfolder, 'DoC_histograms', Name), ...
                    handles.handles.DoCFigPerROI, '-dtiff');

                close gcf
                
                CellData_out(whichPointsInROI, NDatacolumns + 4) = Data_DegColoc1.DoC;
                CellData_out(whichPointsInROI, NDatacolumns + 5) = Data_DegColoc1.Lr;
                CellData_out(whichPointsInROI, NDatacolumns + 6) = Data_DegColoc1.D1_D2;
                CellData_out(whichPointsInROI, NDatacolumns + 7) = Data_DegColoc1.Lr_rAboveThresh;
                CellData_out(whichPointsInROI, NDatacolumns + 8) = Data_DegColoc1.Density;
                
                Data_DoC{roiIter} = Data_DegColoc1;


                % Average density for each region
                DensityROI{roiIter} = [size([CA1;CA2],1)/SizeROI1^2, ...
                    size(CA1,1)/SizeROI1^2, ...
                    size(CA2,1)/SizeROI1^2];

            end
            %AvDensityCell(nt,:)=mean (DensityROI,1);
        end

%     A = Data_DoC(~cellfun('isempty', Data_DoC(:)));
%     DoC1 = cell2mat(cellfun(@(x) x.DoC(x.Ch == 1), A(:), 'UniformOutput', false));
%     DoC2 = cell2mat(cellfun(@(x) x.DoC(x.Ch == 2), A(:), 'UniformOutput', false));
%     DoC1(DoC1==0) = [];
%     DoC2(DoC2==0) = [];
    
    % DoC1 and DoC2 are DoC scores for chan1->2 and chan2->1 for ALL
    % evaluated points above Lr_r threshold
    DoC1CumHist = zeros(100, 1);
    DoC2CumHist = zeros(100, 1);
        
        DoC1CumHist = DoC1CumHist + histc(CellData_out((CellData_out(:,NDatacolumns + 1) > 0) & (CellData_out(:,12) == 1) & (CellData_out(:, NDatacolumns + 7) == 1), NDatacolumns + 4), ...
            linspace(-1, 1, 100));
        DoC2CumHist = DoC2CumHist + histc(CellData_out((CellData_out(:,NDatacolumns + 1) > 0) & (CellData_out(:,12) == 2) & (CellData_out(:, NDatacolumns + 7) == 1), NDatacolumns + 4), ...
            linspace(-1, 1, 100));
  
    DoC1 = DoC1CumHist/sum(DoC1CumHist(:));
    DoC2 = DoC2CumHist/sum(DoC2CumHist(:));

    % Plot summary data
    handles.handles.DoCFig = figure('color', [1 1 1], 'inverthardcopy', 'off');

    handles.handles.DoCAx(1) = subplot(2,1,1);
%     histHand = histogram(handles.handles.DoCAx(1), DoC1, 100);
    histHand = bar(handles.handles.DoCAx(1), linspace(-1, 1, 100), DoC1);
    histHand.FaceColor = Chan1Color;
    histHand.EdgeColor = rgb(52, 73, 94);
    set(handles.handles.DoCAx(1), 'XLim', [-1 1]);
    xlabel(handles.handles.DoCAx(1), 'DoC Score Ch1 -> Ch2', 'Fontsize', 20);
    ylabel(handles.handles.DoCAx(1), 'Frequency','FontSize',20);
    set(handles.handles.DoCAx(1),'FontSize',20)

    handles.handles.DoCAx(2) = subplot(2,1,2);
%     histHand = histogram(handles.handles.DoCAx(2), DoC2, 100);
    histHand = bar(handles.handles.DoCAx(2), linspace(-1, 1, 100), DoC2);
	histHand.FaceColor = Chan2Color;
    histHand.EdgeColor = rgb(52, 73, 94);
    set(handles.handles.DoCAx(2), 'XLim', [-1 1]);
    xlabel(handles.handles.DoCAx(2), 'DoC Score Ch2 -> Ch1', 'Fontsize', 20);
    ylabel(handles.handles.DoCAx(2), 'Frequency','FontSize',20);
    set(handles.handles.DoCAx(2),'FontSize',20)

    drawnow;
  
    % Save the figure
    try 
        print(fullfile(Outputfolder, 'DoC_histograms', 'Pooled_DoC_histogram.tif'), ...
            handles.handles.DoCFig, '-dtiff');
    catch
        currFig = getframe(handles.handles.DoCFig);
        imwrite(currFig.cdata, fullfile(Outputfolder, 'DoC_histograms', 'Pooled_DoC_histogram.tif'));
    end

    close gcf;
      
    save(fullfile(Outputfolder, 'Data_for_Cluster_Analysis.mat'),'Data_DoC','DensityROI'); %ROIData removed!

end
