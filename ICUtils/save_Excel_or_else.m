function save_Excel_or_else(xlsname,captions,datas)
                           
                           if ispc
                                try
                                    if isnumeric(datas) 
                                        datas_ = num2cell(datas); 
                                    elseif iscell(datas)
                                        datas_ = datas;
                                    end
                                    xlswrite( xlsname,[captions; datas_]);
                                catch
                                    disp('can not write output as xls, save as mat file instead');                                    
                                    fname = xlsname(1:length(xlsname)-4);
                                    save([fname '.mat'],'captions','datas');
                                    %
                                    cap = [];
                                    for ci=1:numel(captions)
                                        elem = ['"' char(captions(ci)) '"'];
                                        if ci<numel(captions)
                                            elem = [elem ','];
                                        end
                                        cap = [cap elem];
                                    end
                                    fid = fopen( [fname '.csv'], 'w' );
                                    fprintf( fid, '%s\n', cap);
                                    fclose(fid);
                                    if iscell(datas), datas=cell2mat(datas); end
                                    dlmwrite([fname '.csv'],datas,'-append');
                                    % 
                                end
                            else
                                try
                                    if isnumeric(datas) 
                                        datas_ = num2cell(datas); 
                                    elseif iscell(datas)
                                        datas_ = datas;
                                    end                                    
                                    xlwrite( xlsname,[captions; datas_]);
                                catch
                                    disp('can not write output as xls, save as mat file instead');
                                    fname = xlsname(1:length(xlsname)-4);
                                    save([fname '.mat'],'captions','datas');
                                    %
                                    cap = [];
                                    for ci=1:numel(captions)
                                        elem = ['"' char(captions(ci)) '"'];
                                        if ci<numel(captions)
                                            elem = [elem ','];
                                        end
                                        cap = [cap elem];
                                    end
                                    fid = fopen( [fname '.csv'], 'w' );
                                    fprintf( fid, '%s\n', cap);
                                    fclose(fid);
                                    if iscell(datas), datas=cell2mat(datas); end
                                    dlmwrite([fname '.csv'],datas,'-append');
                                    %                                 
                                end
                           end
end

