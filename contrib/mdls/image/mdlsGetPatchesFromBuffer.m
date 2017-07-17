%
%function X = mdlsGetPatchesFromBuffer(buff_handle,buff_info,indices)
%
%
% Input:
%
%
% 
% Output:
%
%
function X = mdlsGetPatchesFromBuffer(buff_handle,buff_info,indexes)
    if nargin < 2
        error('Usage: buff_handle, buff_info, indices.');
    elseif nargin < 3
        indexes = 0; % all patches
    end
    N = buff_info.patch_count;
    w = buff_info.patch_width;
    M = w^2;
    prec = buff_info.precision;
    bppatch = buff_info.bppatch;
    bpcol = buff_info.bppix * w;
    TOC = buff_info.TOC;
    if indexes == 0 % all patches
        indexes = 1:buff_info.patch_count;
    end
    X = zeros(M,length(indexes));
    NI = length(TOC);
    if buff_info.compact_mode
        for idx=1:length(indexes)
            fidx = indexes(idx);
            % need to determine to which image the index belongs
            j = 1;
            while fidx > TOC(j).patch_offset
                j = j + 1;
                if j > NI
                    break;
                end
            end
            j = j - 1;
            fidx2 = fidx - TOC(j).patch_offset -1 ; 
            % then do the math for the 'condensed index' ...
            end_row_advance = floor(fidx2/TOC(j).patches_per_row)*(w-1);
            fidx2 = bpcol*(fidx - 1 + end_row_advance);
            if fidx2 ~= ftell(buff_handle)
                fseek(buff_handle,fidx2,'bof');
            end
            X(:,idx) = fread(buff_handle, M, prec);
        end        
    else
        % much simpler, much more wasteful
        for i=1:length(indexes)
            fidx=indexes(i);
            fseek(buff_handle,fidx*bppatch,'bof');
            X(:,i) = fread(buff_handle, M, prec);
        end
    end
end