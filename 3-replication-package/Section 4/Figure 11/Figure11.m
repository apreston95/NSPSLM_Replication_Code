clc
clear all
close all

dynare NSPSLM_estimate.mod nolog;

%% Combine plots

% Stacked with opaque black padding
combine_pngs('Figure11A.png','Figure11B.png','Figure11.png','vertical',0,false);


function [outRGB, outA] = combine_pngs(file1, file2, outname, direction, padval, transparent)
% combine_pngs  Combine two PNGs (side-by-side or stacked), with alpha support.
% Usage:
%   [RGB,A] = combine_pngs('a.png','b.png','out.png','horizontal');          % transparent pad
%   [RGB,A] = combine_pngs('a.png','b.png','out.png','vertical',255,false);   % opaque pad (white)
%
% Inputs
%   file1, file2 : paths to input PNGs
%   outname      : output filename (e.g. 'out.png'); pass [] to skip saving
%   direction    : 'horizontal' or 'vertical' (default 'horizontal')
%   padval       : pad colour in [0,255], scalar (default 255 = white)
%   transparent  : logical; true => transparent padding (default true)
%
% Outputs
%   outRGB       : uint8 RGB image
%   outA         : uint8 alpha channel (255=opaque). Empty if transparent=false and no alpha written.

    if nargin < 4 || isempty(direction),   direction   = 'horizontal'; end
    if nargin < 5 || isempty(padval),      padval      = 255; end
    if nargin < 6 || isempty(transparent), transparent = true; end

    [I1, ~, A1] = imread(file1);
    [I2, ~, A2] = imread(file2);

    % Ensure 3-channel RGB
    if size(I1,3) == 1, I1 = repmat(I1, [1 1 3]); end
    if size(I2,3) == 1, I2 = repmat(I2, [1 1 3]); end

    I1 = im2uint8(I1);
    I2 = im2uint8(I2);

    % Normalise/construct alpha as uint8
    if isempty(A1), A1 = uint8(255 * ones(size(I1,1), size(I1,2))); else, A1 = im2uint8(A1); end
    if isempty(A2), A2 = uint8(255 * ones(size(I2,1), size(I2,2))); else, A2 = im2uint8(A2); end

    [h1, w1, ~] = size(I1);
    [h2, w2, ~] = size(I2);

    switch lower(direction)
        case 'horizontal'
            H = max(h1, h2);
            W = w1 + w2;

            outRGB = uint8(ones(H, W, 3) * padval);
            outA   = uint8(255 * ones(H, W));

            % Place left image
            outRGB(1:h1, 1:w1, :) = I1;
            outA(1:h1, 1:w1)      = A1;

            % Place right image
            outRGB(1:h2, w1+1:w1+w2, :) = I2;
            outA(1:h2,  w1+1:w1+w2)     = A2;

        case 'vertical'
            H = h1 + h2;
            W = max(w1, w2);

            outRGB = uint8(ones(H, W, 3) * padval);
            outA   = uint8(255 * ones(H, W));

            % Top image
            outRGB(1:h1, 1:w1, :) = I1;
            outA(1:h1, 1:w1)      = A1;

            % Bottom image
            outRGB(h1+1:h1+h2, 1:w2, :) = I2;
            outA(h1+1:h1+h2, 1:w2)      = A2;

        otherwise
            error('direction must be ''horizontal'' or ''vertical''.');
    end

    % Handle padding transparency
    if transparent
        % Make padding transparent where nothing was written: detect by alpha==255 AND colour==padval
        padMask = all(outRGB == padval, 3);        % padded colour regions
        % If a pixel is from an image, its alpha was set to A1/A2 above; keep it.
        % If it's only padding, set alpha to 0.
        % However, be careful not to blank valid pixels equal to padval by colour:
        % Only zero where padMask AND alpha has not been set by an image write.
        wroteMask = false(H, W);
        % Build wroteMask from alpha writes above: non-padding areas where A came from I1/I2
        % A simple and robust proxy is where outA ~= 255 OR where colour differs from padval
        wroteMask = ~padMask | (outA ~= 255);  %#ok<NASGU> (kept for clarity)
        % Safer: create an occupancy map directly
        occ = false(H, W);
        switch lower(direction)
            case 'horizontal'
                occ(1:h1, 1:w1) = true;
                occ(1:h2, w1+1:w1+w2) = true;
            case 'vertical'
                occ(1:h1, 1:w1) = true;
                occ(h1+1:h1+h2, 1:w2) = true;
        end
        outA(~occ) = 0;  % transparent padding
    else
        % Force fully opaque output (no transparency)
        outA = [];
    end

    if ~isempty(outname)
        if isempty(outA)
            imwrite(outRGB, outname);
        else
            imwrite(outRGB, outname, 'Alpha', outA);
        end
    end
end

