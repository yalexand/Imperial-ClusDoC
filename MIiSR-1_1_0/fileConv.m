function fileConv (inputFile, binSize, Conditions)

%==========================================================================
%                              FUNCTION
% Function convert particle position tables into a matlab-compatible matrix
%
% binSize's of 8000-10,000 are best.
%
%==========================================================================
%                             INPUT/OUTPUT
% Inputs:
%   -inputFile - ASCII file to process
%   -binSize - number of points to process between hard-drive writes.
%              Recommended 8000-10,000.
%   -iFilter - filter based on intensity (# photons)
%   -Conditions - a structure containing:
%       -.iFilter - intensity filter; removes pixels with fewer than
%                   .iFilter photons. Set to 0 to prevent intensity 
%                   filtering.
%       -.cropRange - range (in nm) of image to keep; outside of range is
%                    cropped.  Structured [xmin, xmax, ymin, ymax].  [] = no
%                    cropping
%       -.Scale - pixels -> nm conversion factor.  QuickPALM data is
%                 usually pre-scaled, so this is usually set to 1.
%       -.iScale - value to multiply intensity data in table to convert to
%                  photons
%       -.dataType - type of data table to convert:
%           -1: QuckPALM format
%           -2: Leica GSDR format
%           -3: Tab-delineated X/Y coordinates
%           -4: Tab-delineated X/Y/Intensity coordinates
%           -5: Tab-delineated X/Y/Z coordinates
%           -6: Tab-delineated X/Y/Z/Intensity coordinates
%           -7: Zeiss ELYRA PS1 (*.txt)
%
% Output: a .mat file with the same pre-'.ascii' prefix as the input file, 
%         same data passed to outputMatrix variable if required. Saved with
%         the following columns:
%           -Col 1-3: X/Y/Z coordinates, in nm
%           -Col 4: Precision, in photons (column 2 in quickPALM output)
%           -Col 5: Precision, in nm
%
%==========================================================================
%                       NOTES & DEPENDENCIES
%
%   1) Can be called using only first two options, for a straight
%         conversion without filtering.
%
%   2)  Not parallelized, to allow for use in parfor loops.
%
%   3) Precision is approximated using the formula  p = 200/sqrt(photons),
%     as per www.leica-microsystems.com/fileadmin/downloads/Leica%20SR%20GSD/Application%20Notes/Leica_SR_GSD-ApplicationNote_en.pdf
%
%==========================================================================
%                             CITATION
%
% This script is provided as a supplemental material in:
%
%   Fabiana A. Caetano, Brennan S. Dirk, Joshua H.K. Tam, P. Craig 
%       Cavanagh, Maria Goiko, Stephen S.G. Ferguson, Stephen H. Pasternak,
%       Jimmy D. Dikeakos, John R. de Bruyn, Bryan Heit. MIiSR: Analysis of 
%		Molecular Interactions in Super-Resolution Imaging Enables the Study
%		of Protein Interactions, Dynamics and Formation of Multi-protein 
%		Structures. 2015.  PLoS Computational Biology
% 
% Please reference this paper in any publications which use this script for 
% analysis.
%==========================================================================
%                             REVISIONS
%
% V1.0.1 - remove entries with x/y coordinates of 0

%% Pre-processing
if nargin == 2
    Conditions.Scale = 1;
    Conditions.iFilter = 0;
    Conditions.cropRange = [];
    Conditions.iScale = 1;
    Conditions.dataType = 1;
elseif nargin ~= 3
    error ('fileConv must be called with 2 or 3 input variables.');
end


fID = fopen (inputFile, 'r');

%prepare output file
outputFile = char(inputFile);
pd = max(strfind(outputFile, '.'));
outputFile(pd:end) = [];
outputFile = [outputFile '.mat'];
firstSave = 1;

tMat = fgetl(fID); %ignore first line
endFile = true;
tPos = 0;

%% Load/convert file
display (' ');
display (['Convering ' inputFile '.']);
if Conditions.dataType ~= 2
    while endFile
        fPos = 1;
        tPos = tPos+1;
        %outMat = zeros(binSize,15); %pre-allocate for speed
        while endFile && rem(fPos, binSize+1)
            tMat = fgetl(fID);
            if tMat ~= -1 %end of file
                tMat = rtext(char(tMat), '\t', '', '', 'textsource');
                outMat(fPos,:) = cell2mat(tMat);
                fPos = fPos + 1;
            else
                outMat(fPos:end,:) = []; %remove empty indicies
                endFile = false;
            end
        end %end inner while
        %saveData
        if firstSave
            save  (outputFile, 'outMat', '-v7.3');
            matObj = matfile(outputFile, 'Writable',true);
            firstSave = 0;
        elseif endFile
            matObj.outMat(end+1:end+length(outMat),:) = outMat;
        else %if end of process
            matObj.outMat(end+1:end+length(outMat),:) = outMat;
        end

    end %end outer while
else
    while endFile
        fPos = 1;
        tPos = tPos+1;
        outMat = zeros(binSize,9); %pre-allocate for speed
        while endFile && rem(fPos, binSize+1)
            tMat = fgetl(fID);
            if tMat ~= -1 %end of file
                tMat = rtext(tMat, ',', '', '', 'textsource');
                outMat(fPos,:) = cell2mat(tMat);
                fPos = fPos + 1;
            else
                outMat(fPos:end,:) = []; %remove empty indicies
                endFile = false;
            end
        end %end inner while
        %saveData
        if firstSave
            save  (outputFile, 'outMat', '-v7.3');
            matObj = matfile(outputFile, 'Writable',true);
            firstSave = 0;
        elseif endFile
            matObj.outMat(end+1:end+binSize,1:9) = outMat;
        else %if end of process
            matObj.outMat(end+1:end+length(outMat),1:9) = outMat;
        end

    end %end outer while
    
end
fclose(fID);
clear tMat endFile fID firstSave fPos matObj pd tPos

%% Reformat File
load (outputFile);
if Conditions.dataType == 1
    tMat(:,1:3) = outMat(:,5:7); %x/y/z positions, pre-scaled
    tMat(:,4) = outMat(:,2); %intensity, in # photons
    tMat(:,5) = 0;
elseif Conditions.dataType == 2
    tMat(:,1:2) = outMat(:,4:5); %x/y positions
    tMat(:,3) = outMat(:,7); %z positions
    tMat(:,4) = outMat(:,6); %intensity, in # photons
    tMat(:,5) = 0;
elseif Conditions.dataType == 3
    tMat(:,1:2) = outMat(:,1:2); %x/y positions
    tMat(:,3) = 0; %z positions
    tMat(:,4) = 0; %intensity, in # photons
    tMat(:,5) = 0;
elseif Conditions.dataType == 4
    tMat(:,1:2) = outMat(:,1:2); %x/y positions
    tMat(:,3) = 0; %z positions
    tMat(:,4) = outMat(:,3); %intensity, in # photons
    tMat(:,5) = 0;
elseif Conditions.dataType == 5
    tMat(:,1:3) = outMat(:,1:3); %x/y/z positions
    tMat(:,4) = outMat(:,3); %intensity, in # photons
    tMat(:,5) = 0;
elseif Conditions.dataType == 6 %Conditions.dataType == 6
    tMat(:,1:4) = outMat(:,1:4); %x/y/z/intensity positions
    tMat(:,5) = 0;
else
	tMat(:,1:2) = outMat(:,5:6); %x/y positions
    tMat(:,3) = 0; %z positions
    tMat(:,4) = outMat(:,8); %intensity, in # photons
    tMat(:,5) = outMat(:,7); %precission
end
    
%Scaling 
tMat(:,1:3) = tMat(:,1:3).*Conditions.Scale; %scale x/y/z
if Conditions.dataType ~= 7 %ignore Zeiss data, as is pre-scaled
    tMat(:,4) = tMat(:,4).*Conditions.iScale; %convert intensity to photons
    if max(tMat(:,5)) == 0
        tMat(:,5) = 200./sqrt(tMat(:,4)); %calculate precision
    end
end

%==================== New In V1.0.1 ==================
%check for entries with zero as both x and y coordinates
isZero = find(tMat(:,1) == 0);
if ~isempty(isZero)
    if tMat(isZero,2) == 0
        tMat(isZero,:) = [];
    end
end
%================ End New In V1.0.1 ==================

outMat = tMat;
clear tMat isZero

%intensity filter, if required
if Conditions.iFilter
    outMat = outMat(outMat(:,4) >= Conditions.iFilter,:);
end

%crop image, if required
if ~isempty(Conditions.cropRange)
    outMat = outMat(outMat(:,1) >= Conditions.cropRange(1) & outMat(:,1) <= Conditions.cropRange(2) & outMat(:,2) >= Conditions.cropRange(3) & outMat(:,2) <= Conditions.cropRange(4),:);
end

%% Output & Cleanup
save  (outputFile, 'outMat', '-v7.3');
clear outMat
end % end main function
    
function      [data, result]= rtext(text, delimiter, comment, quotes, options)

% Read input from file, derived from readtext.m, available from: 
%         www.mathworks.com/matlabcentral/fileexchange/10946-readtext

	
	% Read (or set to default) the input arguments.
	if((nargin < 1) || ~ischar(text))		% Is there a file name?
		error('First argument must be the source!'); 
	end
	opts.delimiter=		',';									% Default delimiter value.
	opts.comment=		'';										% Default comment value.
	opts.quotes=		'';										% Default quotes value.
	opts.file=			true;									% Default is file name in source.
	opts.format=		'mixed';								% Default format value.
	opts.op_empty=		[];										% Ignore empties. 
	opts.fname=			'text';
	if(nargin >= 2), opts.delimiter=	delimiter;		end
	if(nargin >= 3), opts.comment=		comment;		end
	if(nargin >= 4), opts.quotes=		quotes;			end
	if(nargin >= 5), options=			options;
	else			 options=			'';
	end
	
	if    (~ischar(opts.delimiter) || isempty(opts.delimiter))
		error('Argument ''delimiter'' must be a non-empty string.');
	elseif(~ischar(opts.comment)   || (length(opts.comment) > 1))
		error('Argument ''comment'' must be a string of maximum one character.');
	elseif(~ischar(opts.quotes)    || (length(opts.quotes) > 2))
		error('Argument ''quotes'' must be a string of maximum two characters.');
	elseif(~ischar(options)   )
		error('Argument ''options'' must be a string.');
	end
	mywaitbar=			@emptywaitbar;							% Default is using no waitbar ...
	th=					[];										% ... so empty waitbar handle.
	options=			lower(options);
	if(~isempty(strfind(options, 'textsource'))),	opts.file= false;			end
	if(~isempty(strfind(options, 'numeric'))),		opts.format= 'numeric';		end
	if(~isempty(strfind(options, 'textual'))),		opts.format= 'textual';		end
	if(~isempty(strfind(options, 'empty2zero'))),	opts.op_empty= 0;	% Replace by zero
	elseif(strcmp(opts.format, 'numeric') || ~isempty(strfind(options, 'empty2nan')))
													opts.op_empty= NaN;	% Replace by NaN.
	end
	if(strcmp(opts.format, 'textual')), opts.op_empty= num2str(opts.op_empty);	end	% Textual 'empty'.
	
	% Set the default return values.
	result.min=		Inf;
	result.max=		0;
	result.quote=	0;
	
	% Read the file.
	if(opts.file)
		opts.fname=		text;
		[fid, errmess]=	fopen(text, 'r');							% Open the file.
		if(fid < 0), error(['Trying to open ''' opts.fname ''': ' errmess]); end
		text=			transpose(fread(fid, 'uchar=>char'));		% Read the file.
		fclose(fid);												% Close the file.
	end
	
	if(~isempty(strfind(options, 'usewaitbar')))
		mywaitbar=		@waitbar;
		th= 			mywaitbar(0, '(readtext) Initialising...');% Show waitbar.
		set(findall(th, '-property', 'Interpreter'), 'Interpreter', 'none');% No (La)TeX formatting. 
	end
	
	% Clean up the text.
	eol=			char(10);									% Using unix-style eol. 
	text=			strrep(text, [char(13) char(10)], eol);		% Replace Windows-style eol.
	text=			strrep(text, char(13), eol);				% Replace MacClassic-style eol.
	if(~isempty(opts.comment))									% Remove comments.
		text=	regexprep(text, ['^\' opts.comment '[^' eol ']*' eol], '');	% Remove entire commented lines. 
		text=	regexprep(text, [ '\' opts.comment '[^' eol ']*'], '');		% Remove commented line endings. 
	end
	if(isempty(text) || text(end) ~= eol),	text= [text eol];	end	% End string with eol, if none.
	
	% Find column and row dividers.
	opts.delimiter=		strrep(opts.delimiter, '\t', char( 9));	% Convert to one char, quicker?
	opts.delimiter=		strrep(opts.delimiter, '\n', char(10));
	opts.delimiter=		strrep(opts.delimiter, '\r', char(13));
	opts.delimiter=		strrep(opts.delimiter, '\f', char(12));
	opts.delimiter=		strrep(opts.delimiter, [char(13) char(10)], eol);% Replace Windows-style eol.
	opts.delimiter=		strrep(opts.delimiter, char(13), eol);	% Replace MacClassic-style eol.
	if(1 == length(opts.delimiter))								% Find column dividers quickly.
		delimS=		find((text == opts.delimiter) | (text == eol));
		delimE=		delimS;
	elseif(isempty(regexp(opts.delimiter, '[\+\*\?\|\[\^\$<>\.\\]', 'once'))) % Find them rather quickly.
		delimS=		strfind(text, opts.delimiter);
		eols=		find(text == eol);
		delimE=		union(eols, delimS + length(opts.delimiter) - 1);
		delimS=		union(eols, delimS);
	else														% Find them with regexp.
		[delimS, delimE]=	regexp(text, [opts.delimiter '|' eol]);
	end
	divRow=			[0, find(text == eol)];						% Find row dividers+last.
	
	% Keep quoted text together.
	if(~isempty(opts.quotes))									% Should we look for quotes?
		if((length(opts.quotes) == 1) || (opts.quotes(1) == opts.quotes(2)))	% Opening char == ending.
			exclE=			find(text == opts.quotes(1));
			exclS=			exclE(1:2:end);
			exclE=			exclE(2:2:end);
		else													% Opening char ~= closing.
			exclS=			find(text == opts.quotes(1));
			exclE=			find(text == opts.quotes(2));
		end
		if((length(exclS) ~= length(exclE)) || any(exclS > exclE))
			close(th);											% Close waitbar or it'll linger.
			error('Opening and closing quotes don''t match in %s.', opts.fname); 
		elseif(~isempty(exclS))									% We do have quoted text.
			mywaitbar(0, th, '(readtext) Doing quotes...');		% Inform user.
			r=		1;
			rEnd=	length(exclS);
			n=		1;
			nEnd=	length(delimS);
			result.quote=	rEnd;
			exclS(end+1)=	0;	% This and next lines are needed in cases with only one quote in text.
			exclE(end+1)=	0;
			while((n <= nEnd) && (r <= rEnd)) % "Remove" delimiters and newlines within quotes.
				while((r <= rEnd) && (delimS(n) > exclE(r))), r= r+1;	end	% Next end-quote after newline.
				while((n <= nEnd) && (delimS(n) < exclS(r))), n= n+1;	end	% Next newline after strart-quote.
				while((n <= nEnd) && (delimS(n) >= exclS(r)) && (delimS(n) <= exclE(r)))
					delimS(n)=	0;											% Newlines inside quote. 
					n=			n+1;
				end
				mywaitbar(n/nEnd);								% Update waitbar.
			end
			mywaitbar(1);
			delimE=	delimE(delimS > 0);
			delimS=	delimS(delimS > 0);
		end
	end
	delimS=		delimS-1;										% Last char before delimiter.
	delimE=		[1 delimE(1:end-1)+1];							% First char after delimiter.
	
	% Do the stuff: convert text to cell (and maybe numeric) array.
	mywaitbar(0, th, sprintf('(readtext) Processing ''%s''...', opts.fname));
	r=				1;
	c=				1;											% Presize data to optimise speed.
	data=			cell(length(divRow), ceil(length(delimS)/(length(divRow)-1)));
	nums=			zeros(size(data));							% Presize nums to optimise speed.
	nEnd=			length(delimS);								% Prepare for a waitbar.
	istextual=		strcmp(opts.format, 'textual');
	for n=1:nEnd
		temp= 			strtrim(text(delimE(n):delimS(n)));		% Textual item.
		data{r, c}= 	temp;									% Textual item.
		if(~istextual)
			lenT=			length(temp);
			if(lenT > 0)
				% Try to get 123, 123i, 123i + 45, or 123i - 45
				[a, count, errmsg, nextindex]=	sscanf(temp, '%f %1[ij] %1[+-] %f', 4);
				if(isempty(errmsg) && (nextindex > lenT))
					if    (count == 1),			nums(r, c)=		a;
					elseif(count == 2),			nums(r, c)=		a(1)*i;
					elseif(count == 4),			nums(r, c)=		a(1)*i + (44 - a(3))*a(4);
					else						nums(r, c)=		NaN;
					end
				elseif(regexpi(temp, '[^0-9eij \.\+\-]', 'once'))	% Remove non-numbers.
					nums(r, c)=		NaN;
				else
					nums(r, c)= 	s2n(temp, lenT);			% Some other kind of complex number.
				end
			else			nums(r, c)=	NaN;
			end
		end
		if(text(delimS(n)+1) == eol)							% Next row.
			result.min=		min(result.min, c);					% Find shortest row.
			result.max=		max(result.max, c);					% Find longest row.
			r=				r+1;
			c=				0;
			if(bitand(r, 15) == 0), mywaitbar(n/nEnd);	end		% Update waitbar.
		end
		c=				c+1;
	end
	
	% Clean up the conversion and do the result statistics.
	mywaitbar(0, th, '(readtext) Cleaning up...');				% Inform user.
	data=				data(1:(r-1), 1:result.max);			% In case we started off to big.
	if(~strcmp(opts.format, 'textual')), nums= nums(1:(r-1), 1:result.max);	end		% In case we started off to big.
	if(all(cellfun('isempty', data)))
		data= 			{};
		r=				1;
		result.min=		0;
		result.max=		0;
	else
		while((size(data, 2) > 1) && all(cellfun('isempty', data(end, :))))	% Remove empty last lines. 
			data=	data(1:end-1, :); 
			nums=	nums(1:end-1, :); 
			r=		r-1;
		end 
		while((size(data, 1) > 1) && all(cellfun('isempty', data(:, end))))	% Remove empty last columns. 
			data=	data(:, 1:end-1); 
			nums=	nums(:, 1:end-1); 
			c=		c-1;
		end 
	end
	result.rows=		r-1;
	empties=			cellfun('isempty', data);				% Find empty items.
	result.emptyMask=	empties;
	if(strcmp(opts.format, 'textual'))
		result.numberMask=	repmat(false, size(data));			% No numbers, all strings.
		result.stringMask=	~empties;							% No numbers, all strings.
		data(empties)=		{opts.op_empty};					% Set correct empty value.
	else
		if(isempty(data)),	result.numberMask=	[];
		else				result.numberMask=	~(isnan(nums) & ~strcmp(data, 'NaN'));	% What converted well.
		end
		if(strcmp(opts.format, 'numeric'))
			nums(empties)=		opts.op_empty;					% Set correct empty value.
			data=				nums;							% Return the numeric array.
			result.stringMask=	~(empties | result.numberMask);	% Didn't convert well: so strs.
		else
			data(result.numberMask)= num2cell(nums(result.numberMask));	% Copy back numerics.
			data(empties)=		{opts.op_empty};						% Set correct empty value.
			result.stringMask=	cellfun('isclass', data, 'char');	% Well, the strings.
		end
	end
	result.empty=		sum(result.emptyMask(:));				% Count empties.
	result.numberMask=	result.numberMask & ~result.emptyMask;	% Empties don't count.
	result.number=		sum(result.numberMask(:));				% Count numbers.
	result.stringMask=	result.stringMask & ~result.emptyMask;	% Empties don't count.
	result.string=		sum(result.stringMask(:));				% Count strings.
	
 	close(th);													% Removing the waitbar. 
end


function x= s2n(s, lenS)
	x=		NaN;

	% Try to get 123 + 23i or 123 - 23i
	[a,count,errmsg,nextindex] = sscanf(s,'%f %1[+-] %f %1[ij]',4);
	if(isempty(errmsg) && (nextindex > lenS))
		if(count == 4),				x=		a(1) + (44 - a(2))*a(3)*i;
		end
		return
	end

	% Try to get i, i + 45, or i - 45
	[a,count,errmsg,nextindex] = sscanf(s,'%1[ij] %1[+-] %f',3);
	if(isempty(errmsg) && (nextindex > lenS))
		if(count == 1),				x=		i;
		elseif(count == 3),			x=		i + (44 - a(2))*a(3);
		end
		return
	end

	% Try to get 123 + i or 123 - i
	[a,count,errmsg,nextindex] = sscanf(s,'%f %1[+-] %1[ij]',3);
	if(isempty(errmsg) && (nextindex > lenS))
		if(count == 1),				x=		a(1);
		elseif(count == 3),			x=		a(1) + (44 - a(2))*i;
		end
		return
	end

	% Try to get -i, -i + 45, or -i - 45
	[a,count,errmsg,nextindex] = sscanf(s,'%1[+-] %1[ij] %1[+-] %f',4);
	if(isempty(errmsg) && (nextindex > lenS))
		if(count == 2),				x=		(44 - a(1))*i;
		elseif(count == 4),			x=		(44 - a(1))*i + (44 - a(3))*a(4);
		end
		return
	end

	% Try to get 123 + 23*i or 123 - 23*i
	[a,count,errmsg,nextindex] = sscanf(s,'%f %1[+-] %f %1[*] %1[ij]',5);
	if(isempty(errmsg) && (nextindex > lenS))
		if(count == 5),				x=		a(1) + (44 - a(2))*a(3)*i;
		end
		return
	end

	% Try to get 123*i, 123*i + 45, or 123*i - 45
	[a,count,errmsg,nextindex] = sscanf(s,'%f %1[*] %1[ij] %1[+-] %f',5);
	if(isempty(errmsg) && (nextindex > lenS))
		if(count == 1),				x=		a;
		elseif(count == 3),			x=		a(1)*i;
		elseif(count == 5),			x=		a(1)*i + (44 - a(4))*a(5);
		end
		return
	end

	% Try to get i*123 + 45 or i*123 - 45
	[a,count,errmsg,nextindex] = sscanf(s,'%1[ij] %1[*] %f %1[+-] %f',5);
	if(isempty(errmsg) && (nextindex > lenS))
		if(count == 1),				x=		i;
		elseif(count == 3),			x=		i*a(3);
		elseif(count == 5),			x=		i*a(3) + (44 - a(4))*a(5);
		end
		return
	end

	% Try to get -i*123 + 45 or -i*123 - 45
	[a,count,errmsg,nextindex] = sscanf(s,'%1[+-] %1[ij] %1[*] %f %1[+-] %f',6);
	if(isempty(errmsg) && (nextindex > lenS))
		if(count == 2),				x=		(44 - a(1))*i;
		elseif(count == 4),			x=		(44 - a(1))*i*a(4);
		elseif(count == 6),			x=		(44 - a(1))*i*a(4) + (44 - a(5))*a(6);
		end
		return
	end

	% Try to get 123 + i*45 or 123 - i*45
	[a,count,errmsg,nextindex] = sscanf(s,'%f %1[+-] %1[ij] %1[*] %f',5);
	if(isempty(errmsg) && (nextindex > lenS))
		if(count == 5),				x=		a(1) + (44 - a(2))*i*a(5);
		end
		return
	end

	% None of the above cases.
	x=		NaN;

end


function emptywaitbar(varargin)
end
