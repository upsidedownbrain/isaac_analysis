function [ tseries, tseries_shape, rois ] = isaac_adapt_files( varargin )
% Helper function to adapt nifti files for ISAAC analysis. It can be used
% to adapt the time series, the regions of interest or both. For it to
% work, load_nifti() funciton from the freesurfer matlab toolbox is
% necessary
%
% Example:
%
%   
%   file_ts = 'bold.nii.gz'
%
%   % we might have 20 rois, and want to compare the first 10 with the
%   % latter 10
%   files_rois = {'roi01.nii.gz', 'roi02.nii.gz', ..... 'roi20.nii.gz'};
%   roi_names = {'region01', 'region02', .... 'region20'};
%   roi_pairs = struct;
%   roi_pairs.as_x = 1:10;
%   roi_pairs.as_y = 11:20;
%
%   [ ts, sz, rois ] = isaac_adapt_files( 'tseries_file', file_ts,...
%       'roi_files', files_rois, 'roi_names', roi_names, 'roi_pairs', roi_pairs );
% 
% It gets the file name of the fmri nifti and returns the time series as a
% NTimePoints x NVoxels array.



% --------------- parse args ----------------------------------------------
if rem(length(varargin), 2)
    warning('Odd number of optional arguments. Optional arguments must be in name-value pairs');
end
% Set default values:
tseries_file = [];
roi_files = [];
roi_names = [];
roi_pairs = [];

% Get user defined params:
for k = 1:2:length(varargin)-1
    name = varargin{k};
    switch name
        case 'tseries_file';
            tseries_file = varargin{k+1};
        case 'roi_files';
            roi_files = varargin{k+1};
        case 'roi_names';
            roi_names = varargin{k+1};
        case 'roi_pairs';
            roi_pairs = varargin{k+1};
        otherwise
            warning('%s: Argument name %s not known, it will be ignored', ...
                mfilename, name);
    end
end

% detect the file type (nifti or cifti)
if isempty(tseries_file) && isempty(roi_files)
    error('No files to adapt were specified.');
end



% ------- load and adapt tseries ------------------------------------------
tseries = [];
tseries_shape = [];
rois = [];

switch get_file_type(tseries_file)
    case 'nifti'
        tseries_format = 'nifti';
        [tseries, tseries_shape] = nii_adapt_tseries(tseries_file);
    case 'cifti'
        error('%s: cifti class not implemented yet', mfilename);
        tseries_format = 'cifti';
        [tseries, tseries_shape] = cifti_adapt_tseries(tseries_file);
    otherwise
        error('%s: unsupported file format for input time series', mfilename);
end


% ------- load and adapt rois ---------------------------------------------
switch get_file_type(roi_files)
    case 'nifti'
        rois = nii_adapt_rois(roi_files, roi_names, roi_pairs);
        rois_format = 'nifti';
    case 'cifti'
        error('cifti class not implemented yet, try reading and addapt the files yourself');
        rois_format = 'cifti';
        rois = cifti_adapt_rois(roi_files, roi_names, roi_pairs);
    otherwise
        error('Unsupported file format for input rois, must be ALL nifti files.');
end



% if both rois and tseries were specified, check they are compatible
if ~isempty(roi_files) && ~isempty(tseries_file)
    % check they are both of the same class:
    if ~strcmp(rois_format, tseries_format)
        error('tseries and roi files are not of the same type. They must be both in nifti format or both in cifti format');
    end
    
    % check the shapes
    for k = 1:length(rois.shape)
        if any(rois.shape{k} - tseries_shape(1:3) ~= 0)
            error('ROI isn''t in the same space as time series. \n\troi file:%s\n\n\troi shape:(%d,%d,%d)\n\n\ttseries shape:(%d,%d,%d,%d)\n', ...
                roi_files{k}, rois.shape{k}, tseries_shape);
        end
    end
end




function file_type = get_file_type(file_in)

is_cifti = regexp(file_in, '\.(p|d|pd|dp)((scalar)|(tseries)|(label))\.nii$');
is_nifti = regexp(file_in, '.nii(\.gz){0,1}$');

if iscellstr(file_in)
    is_cifti = ~cellfun(@isempty, is_cifti);
    is_nifti = ~cellfun(@isempty, is_nifti) & ~is_cifti;
elseif ischar(file_in)
    is_cifti = ~isempty(is_cifti);
    is_nifti = ~isempty(is_nifti) & ~is_cifti;
else
    file_type = '';
    return;
end

if all(is_nifti)
    file_type = 'nifti';
elseif all(is_cifti)
    file_type = 'cifti';
else
    file_type = '';
end


% -------- actual addapting functions -------------------------------------
function [tseries, tseries_shape] = nii_adapt_tseries(tseries_file)
% This function will assume that data is a mxnxpxT matrix, where T is the
% time dimension.
% 
% This will also cast the volume to double class.

% read the file:
nii = load_nifti(tseries_file);

% Reshape into a #voxels x # #volumes matrix:
tseries_shape = size(nii.vol);
tseries = permute(reshape(nii.vol, prod(tseries_shape(1:3)), tseries_shape(4)), [2,1]);


function [rois, tseries_shape] = nii_adapt_rois(roi_files, roi_names,  roi_pairs)


% if the names are not specified, get them from the file names:
if (numel(roi_files)>1 && numel(roi_names)~=numel(roi_files)) || isempty(roi_names)
    warning('The roi names will be specified automatically.');
    roi_names = cellfun(@(c) strsplit(c, filesep), roi_files, 'uni', false);
    roi_names = cellfun(@(c) c{end}, roi_names, 'uni', false);
end

% read the file:
rois = struct;
rois.names = roi_names;
rois.idx   = cell(size(roi_files));
rois.shape = cell(size(roi_files));

if ~isempty(roi_pairs)
    rois.as_x  = roi_pairs.as_x;
    rois.as_y  = roi_pairs.as_y;
else
    rois.as_x  = 1:numel(roi_files);
    rois.as_y  = 1:numel(roi_files);
end

% this is the actual adaptation of the volume to idx. It is assumed that
% rois are 3d-vols with 0s and 1s.
for k = 1:numel(roi_files)
    nii = load_nifti(roi_files{k});
    
    rois.shape{k} = size(nii.vol);
    rois.idx{k} = find(nii.vol);
end

% store also the file names in the final structure, maybe for debugging
rois.filenames = roi_files;





