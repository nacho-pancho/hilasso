%
% given a list of images and a patch width, save the patches as a long
% continuous file where the patches can be accessed efficiently.
%
% This version saves a B/W version of the patches.
%
% The routine can work in two modes: normal or compact.
% In normal mode, every possible patch of the specified width is stored
% sequentially in the file. This explodes the number of pixels in the original
% image by width^2 so be careful. This is specially useful when working on
% single (or a few) images where one wants maximum access speed and 
% to be able to modify the patches and re-write them online.
%
% The other mode, suitable for building a large database, saves space by using the
% fact that patches are stored column-wise, so to reconstruct any patch in a row
% it is enough to store the non-overlapping contiguous patches in order, so that
% an intermediate overlapping patch can be reconstructed using the last few columns
% of one patch, and the first few of the next. Since the file is a plain byte
% file, this reduces to computing the correct offset each time.
%
% Input:
% 
% image_dir ...... directory where the images are stored.
% image_list ..... which images are included in the cache
% patch_width .... width of the patches
% buffer_prefix .. file name prefix for the generated files. These are
%                  <buffer_prefix>-patches.bin ... binary file with the patch data
%                  <buffer_prefix>-info.mat ....... mat-format info file.
%                  <buffer_prefix>-info.txt ....... text-format info file.
%                  <buffer_prefix>-dc.bin ........ optional, containing the DC of the patches.
% remove_dc ...... If true, removes the DC from each patch and places it in the
%                  same position in the file <buffer_prefix>-dc.bin. Defaults to false.
% compact_mode ... See above. Defaults to false.
%
% Output:
% 
% buffer_info .... Struct with information about the buffer.
%    buffer_info.TOC ...........
%    buffer_info.buffer_prefix . prefix for the name of the buffer and related files
%    buffer_info.byte_count .... total number of bytes written
%    buffer_info.patch_count ... total number of patches written
%    buffer_info.remove_dc ..... wether DC was removed from patches
%    buffer_info.compact_mode .. wether buffer was crated in compact mode
%    buffer_info.bppatch ....... bytes pero 
%    buffer_info.bppix ......... bytes per pixel
%    buffer_info.precision ..... precision with which pixels were stored in Matlab naming convention
% 
function buffer_info = mdlsCreatePatchBuffer(image_dir, ...
                                     image_list, ...
                                     patch_width, ...
                                     buffer_prefix, ...
                                     remove_dc, ...
                                     compact_mode)
    if nargin < 4
        error('Usage: dir, list, width, buffer_name.');
    end
    if ~exist('remove_dc','var')
        remove_dc = false;
    end
    if ~exist('compact_mode','var')
        compact_mode = false;        
    end
    % cannot use compact mode if DC is removed
    if remove_dc
        compact_mode = false;
    end
    NI = length(image_list);
    % TOC
    TOC = struct();
    %
    % init
    %
    if isequal(lower(buffer_prefix(end-3:end)),'.mat')
        warning('Buffer prefix must not contain .mat at the end. Removing.');
        buffer_prefix = buffer_prefix(1:end-4);
    end
    buffer_prefix = sprintf('%s_w%02d',buffer_prefix,patch_width);
    buffer_name = [buffer_prefix '-patches.bin'];
    fh = fopen(buffer_name,'w');
    [em,en]=ferror(fh);
    if (en ~= 0) 
        error(['Error opening ' buffer_name ' for writing: ' em]);
    end
    if remove_dc
        precision = 'int16';
        dc_name = [buffer_prefix '-dc.bin'];
        fdc = fopen(dc_name,'w');
        [em,en]=ferror(fdc);
        if (en ~= 0) 
            error(['Error opening ' dc_name ' for writing: ' em]);
        end
    else
        precision = 'uint8';
    end
    
    byte_count = 0;
    patch_count = 0;
    w = patch_width;
    M = patch_width^2;
    overlap = w - 1;
    if remove_dc
        bppix = 2;
    else
        bppix = 1;
    end
    bppatch = bppix*M;
    %
    % main loop
    %
    finfo = fopen([buffer_prefix '-info.txt'],'w');
    fprintf(finfo,'# NUMIMG\tCOMPACT\tREMDC\tPWIDTH\tBPATCH\n');
    fprintf(finfo,'%8d\t%1d\t%1d\t%2d\t%4d\n',NI,compact_mode,remove_dc,w,bppatch);
    fprintf(finfo,'# ROWS\tCOLS\tOFFSET\t\tBYTES\t\tB.P.ROW\t');
    fprintf(finfo,'PATCHES\t\tPATCHOFFSET\tP.P.ROW\tIMGNAME\t\tIMGPATH\n');    
    for i = 1:NI
        fprintf('Processing %05d/%05d ...\n',i,NI);
        image_name = [image_dir '/' image_list{i} ];
        I=imread(image_name);
        [m,n]=size(I);
        if (size(I,3) == 3) 
            I=rgb2gray(I);
        end
        num_patches = (m-w+1)*(n-w+1);
        if compact_mode
            bytes_written = 0;
            for j=1:(m-w)+1
                bytes_written = bytes_written + bppix * fwrite(fh,I(j:(j+w-1),:),precision);            
            end
        else
            if ~remove_dc
                X      = mdlsDeconstructFast(I,w,overlap);
            else
                [X,DC] = mdlsDeconstruct(I,w,overlap);
            end
            bytes_written = bppix * fwrite(fh,X,precision);
        end
        [em,en]=ferror(fh);
        if (en ~= 0) 
            error(['Error writing to ' buffer_prefix '-patches.bin: ' em]);
        end
        if remove_dc
            fwrite(fdc,DC,'single');
            [em,en]=ferror(fdc);
            if (en ~= 0) 
                error(['Error writing to ' buffer_prefix '-dc.bin: ' em]);
            end
        end
        %
        % record for MAT-format TOC
        %
        TOC(i).m = m;
        TOC(i).n = n;
        TOC(i).offset = byte_count;
        TOC(i).patch_offset = patch_count;
        TOC(i).image_name = image_list{i};
        TOC(i).image_path = image_name;
        TOC(i).num_patches = num_patches;
        TOC(i).num_bytes = bytes_written;
        TOC(i).patches_per_row = TOC(i).n - w + 1;
        TOC(i).bytes_per_row = TOC(i).patches_per_row * M;
        %
        % record for TXT-format TOC
        %
        fprintf(finfo,'%d\t%d\t',TOC(i).m,TOC(i).n);
        fprintf(finfo,'%10d\t%10d\t%4d\t',TOC(i).offset,TOC(i).num_bytes,TOC(i).bytes_per_row);
        fprintf(finfo,'%9d\t%9d\t%4d\t',TOC(i).patch_offset,TOC(i).num_patches,TOC(i).patches_per_row);
        fprintf(finfo,'%12s\t%s\n',TOC(i).image_name,TOC(i).image_path);
        byte_count = byte_count + bytes_written;
        patch_count = patch_count + num_patches;
    end
    %
    % finish
    %
    fclose(finfo);
    fclose(fh);
    if remove_dc
        fclose(fdc);
    end
    buffer_info = struct();
    buffer_info.TOC = TOC;
    buffer_info.patch_width = patch_width;
    buffer_info.patch_dim = M;
    buffer_info.buffer_prefix = buffer_prefix;
    buffer_info.byte_count = byte_count;
    buffer_info.patch_count = patch_count;
    buffer_info.remove_dc = remove_dc;
    buffer_info.compact_mode = compact_mode;
    buffer_info.bppatch = bppatch;
    buffer_info.bppix = bppix;
    buffer_info.precision = precision;

    info_file = [ buffer_prefix '-info.mat' ];
    save(info_file, 'buffer_info');
end