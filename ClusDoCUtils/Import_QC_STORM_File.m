
function [Data,Head_text,Body_text] = Import_QC_STORM_File(fname,handles)

% parameters are 
% peak intensity (photon), 
% x (pixel), 
% y (pixel), 
% z (nm), 
% PSFSigmaX (pixel), 
% PSFSigmaY (pixel), 
% Total intensity (photon), 
% background (photon), 
% SNR (peak to background e-), 
% CRLBx (nm), 
% CRLBy (nm), 
% frame

ParaNum = 12; 

fid=fopen(fname,'rb');
loc=fread(fid,inf,'float');
fclose(fid);

Len=floor(length(loc)/ParaNum);
LocArry=zeros(ParaNum,Len);

LocArry(:)=loc(1:ParaNum*Len);
LocArry=LocArry';

pos=LocArry(:,1)~=0;
LocArry=LocArry(pos,:);

LocArry=sortrows(LocArry,ParaNum);
% figure;
% plot(LocArry(:,2)+0.5,LocArry(:,3)+0.5,'x');
% grid on;

            % pretend that it is WindSTORM ...
            Nloc = length(squeeze(LocArry(:,2)));
            Body_text = cell(1,5);
            Body_text{1} = (1:Nloc)';
            Body_text{2} = LocArry(:,12);
            Body_text{3} = LocArry(:,2);
            Body_text{4} = LocArry(:,3);
            Body_text{5} = LocArry(:,7); % total intensity

            Head_text = '"id","frame","x [pix]","y [pix]","intensity [photon]"';        

            pixelSizenm = handles.pixelSizenm;
            psf_width = handles.WindSTORM_Sigmapix*handles.pixelSizenm;

                Data = [(1:length(Body_text{1}))', ... % index
                        Body_text{2}, ... % frame                    
                        ones(length(Body_text{1}), 1), ...
                        zeros(length(Body_text{1}), 1), ...
                        pixelSizenm*Body_text{3}, ...
                        pixelSizenm*(max(Body_text{4})-Body_text{4}), ...
                        zeros(length(Body_text{1}), 1), ... % precision (uncertainty_xy)
                        Body_text{5}, ... % intensity
                        zeros(length(Body_text{1}), 1), ... % bckgstd
                        ones(length(Body_text{1}), 1), ... % chi square
                        psf_width*ones(length(Body_text{1}), 1), ... % psf width [nm] MAY BE DONE BETTER IF NEEDED
                        ones(length(Body_text{1}), 1), ...
                        ones(length(Body_text{1}), 1)];                    
end
