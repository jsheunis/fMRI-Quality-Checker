% function output = calculateQCmeasures(functional4D_fn, structural_fn, fwhm, spm_dir)
% Function to calculate multiple quality control measures for an fMRI time
% series.
% Steps include preprocessing the raw fMRI and structural data
% (coregistering structural image to first functional image, segmenting the
% coregistered structural image into tissue types, and reslicing the
% segments to the functional resolution image grid), and calling functions
% for standard set of QC measures based on various sources (including
% mainly PC-QAP. Makes use of spm12 batch routines.
% If spm12 batch parameters are not explicitly set, defaults are assumed.
%
% INPUT:
% funcional4D_fn     - filename of pre-real-time functional scan
% structural_fn      - filename of T1-weighted structural scan
% fwhm               - kernel size for smoothing operations
%
% OUTPUT:
% output            - structure with filenames and data
%__________________________________________________________________________
% Copyright (C) Stephan Heunis 2018


% User defined variables
% -------------------------------------------------------------------------
% data_dir = '/Users/jheunis/Documents/MATLAB/rt_fMRI/Data/Rolf finger tapping'; % e.g. '/users/me/matlab/data/subj1'
% functional4D_fn = [data_dir filesep 'real_time_fMRI_WIP_FE_EPI_tapping_SENSE_12_1.nii']; % e.g. [data_dir filesep 'rest.nii']
% structural_fn = [data_dir filesep 'real_time_fMRI_WIP_T1W_3D_TFE_SENSE_7_1.nii']; % e.g. [data_dir filesep 'mprage.nii']
subject = '0051210';
data_dir = '/Users/jheunis/Documents/MATLAB/rtqc_jsh/rtqc_data/0051210'; % e.g. '/users/me/matlab/data/subj1'
functional4D_fn = [data_dir filesep 'rest_1/rest.nii']; % e.g. [data_dir filesep 'rest.nii']
structural_fn = [data_dir filesep 'anat_1/mprage.nii']; % e.g. [data_dir filesep 'mprage.nii']
spm_dir = '/Users/jheunis/Documents/MATLAB/spm12'; % e.g. '/users/me/matlab/spm12'
fwhm = 6; % for preprocessing smoothing steps
intensity_scale = [-6 6]; % scaling for plot image intensity, see what works
% -------------------------------------------------------------------------

% Get image information
func_spm = spm_vol(functional4D_fn);
tsize = size(func_spm);
Nt = tsize(1);
Ni= func_spm(1).dim(1);
Nj= func_spm(1).dim(2);
Nk= func_spm(1).dim(3);

% Preprocess structural and functional data.
[d, f, e] = fileparts(structural_fn);
if exist([d filesep 'rc1' f e], 'file')
    % If a resliced grey matter image is present in the specified data
    % directory, assume that all preprocessing has been completed by
    % standardPreproc, and declare appropriate variables.
    disp('Preproc done!')
    preproc_data = struct;
    preproc_data.forward_transformation = [d filesep 'y_' f e];
    preproc_data.inverse_transformation = [d filesep 'iy_' f e];
    preproc_data.gm_fn = [d filesep 'c1' f e];
    preproc_data.wm_fn = [d filesep 'c2' f e];
    preproc_data.csf_fn = [d filesep 'c3' f e];
    preproc_data.bone_fn = [d filesep 'c4' f e];
    preproc_data.soft_fn = [d filesep 'c5' f e];
    preproc_data.air_fn = [d filesep 'c6' f e];
    preproc_data.rstructural_fn = [d filesep 'r' f e];
    preproc_data.rgm_fn = [d filesep 'rc1' f e];
    preproc_data.rwm_fn = [d filesep 'rc2' f e];
    preproc_data.rcsf_fn = [d filesep 'rc3' f e];
    preproc_data.rbone_fn = [d filesep 'rc4' f e];
    preproc_data.rsoft_fn = [d filesep 'rc5' f e];
    preproc_data.rair_fn = [d filesep 'rc6' f e];
    
    [d, f, e] = fileparts(functional4D_fn);
    preproc_data.rfunctional_fn = [d filesep 'r' f e];
    preproc_data.srfunctional_fn = [d filesep 'sr' f e];
    preproc_data.sfunctional_fn = [d filesep 's' f e];
    preproc_data.mp_fn = [d filesep 'rp_' f '.txt'];
    preproc_data.MP = load(preproc_data.mp_fn);
else
    % If a resliced grey matter image is NOT present in the specified data
    % directory, run standardPreproc
    preproc_data = standardPreproc(functional4D_fn, structural_fn, fwhm, spm_dir);
end

% Calculate brain mask matrices for GM, WM, CSF, and all together
mask_threshold = 0.1;
[GM_img_bin, WM_img_bin, CSF_img_bin] = createBinarySegments(preproc_data.rgm_fn, preproc_data.rwm_fn, preproc_data.rcsf_fn, mask_threshold);
I_GM = find(GM_img_bin);
I_WM = find(WM_img_bin);
I_CSF = find(CSF_img_bin);
mask_reshaped = GM_img_bin | WM_img_bin | CSF_img_bin;
I_mask = find(mask_reshaped);
Nmaskvox = numel(I_mask);

% Detrend 4D time series
F4D_detrended = detrend4D(preproc_data.sfunctional_fn);
F2D_detrended = reshape(F4D_detrended, Ni*Nj*Nk, Nt);

% Statistical measures
F2D_mean = mean(F2D_detrended, 2);
F2D_stddev = std(F2D_detrended, 0, 2);
F2D_zstat = (F2D_detrended - F2D_mean)./F2D_stddev;
F2D_zstat(isnan(F2D_zstat))=0;
F2D_zstat_mean = mean(abs(F2D_zstat),1);
Zstat_mean = mean(F2D_zstat_mean);
F2D_var = var(F2D_detrended,0,2);
F2D_psc = 100*(F2D_detrended./repmat(F2D_mean, 1, Nt)) - 100;
F2D_psc(isnan(F2D_psc))=0;

% Framewise displacement
r = 50; % mm
FD_threshold = 1; % mm
FD_measures = calculateFD(preproc_data.MP, r, FD_threshold);

% DVARS
F2D_diff = [zeros(1, Ni*Nj*Nk); diff(F2D_detrended')]';
DVARS = var(F2D_diff);

% The plot (with FD, DVARS, and mean Zscore per volume)
GM_img = F2D_psc(I_GM, :);
WM_img = F2D_psc(I_WM, :);
CSF_img = F2D_psc(I_CSF, :);
all_img = [GM_img; WM_img; CSF_img];
line1_pos = numel(I_GM);
line2_pos = numel(I_GM) + numel(I_WM);
tf = figure;
fontsizeL = 14;
fontsizeM = 11;
ax1 = subplot(7,1,4:7);
imagesc(ax1, all_img); colormap(gray); caxis(intensity_scale);
title(ax1, 'thePlotSpm','fontsize',fontsizeL)
ylabel(ax1, 'Voxels','fontsize',fontsizeM)
xlabel(ax1, 'fMRI volumes','fontsize',fontsizeM)
hold on; line([1 Nt],[line1_pos line1_pos],  'Color', 'b', 'LineWidth', 2 )
line([1 Nt],[line2_pos line2_pos],  'Color', 'r', 'LineWidth', 2 )
hold off;
ax2 = subplot(7,1,1);
plot(ax2, FD_measures.FD, 'LineWidth', 2); grid;
set(ax2,'Xticklabel',[]);
title(ax2, 'FD','fontsize',fontsizeL)
ylabel(ax2, 'mm','fontsize',fontsizeM)
ax3 = subplot(7,1,2);
plot(ax3, DVARS, 'LineWidth', 2); grid;
set(ax3,'Xticklabel',[]);
title(ax3, 'DVARS','fontsize',fontsizeL)
ylabel(ax3, 'a.u.','fontsize',fontsizeM)
ax4 = subplot(7,1,3);
plot(ax4, F2D_zstat_mean, 'LineWidth', 2); grid;
set(ax4,'Xticklabel',[]);
title(ax4, 'Z-score','fontsize',fontsizeL)
ylabel(ax4, 'a.u.','fontsize',fontsizeM)
print(tf, 'timeseries_summary', '-dpng')

% tSNR
tSNR_2D = F2D_mean./F2D_stddev;
tSNR_brain = mean(tSNR_2D(I_mask));
tSNR_GM = mean(tSNR_2D(I_GM));
tSNR_WM = mean(tSNR_2D(I_WM));
tSNR_CSF = mean(tSNR_2D(I_CSF));

% Metrics
disp(['Number of volumes classified as outliers based on FD>=' num2str(FD_threshold) 'mm: ' num2str(numel(FD_measures.FD_outliers_ind))])
disp(['Total FD: ' num2str(FD_measures.FD_sum)])
disp(['Mean FD: ' num2str(FD_measures.FD_mean)])
disp(['Mean Zscore: ' num2str(Zstat_mean)])
disp(['tSNR (brain): ' num2str(tSNR_brain)])
disp(['tSNR (GM): ' num2str(tSNR_GM)])
disp(['tSNR (WM): ' num2str(tSNR_WM)])
disp(['tSNR (CSF): ' num2str(tSNR_CSF)])

% 3D and 4D images
mask_3D = reshape(mask_reshaped, Ni, Nj, Nk);
tSNR_3D = reshape(tSNR_2D, Ni, Nj, Nk);
F3D_mean = reshape(F2D_mean, Ni, Nj, Nk);
F3D_var = reshape(F2D_var, Ni, Nj, Nk);
F3D_stddev = reshape(F2D_stddev, Ni, Nj, Nk);

% montage3 = createMontage(F3D_var, 5, 1, 'Variance (whole image)');


tSNR_2D_masked = zeros(Ni*Nj*Nk, 1);
tSNR_2D_masked(I_mask, :) = tSNR_2D(I_mask, :);
tSNR_3D_masked = reshape(tSNR_2D_masked, Ni, Nj, Nk);

montage2 = createMontage(F3D_mean, 5, 1, 'Mean EPI (whole image)', 'gray');
print(montage2.f, 'mean_epi', '-dpng')
montage3 = createMontage(F3D_stddev, 5, 1, 'Standard deviation (whole image)', 'parula');
print(montage3.f, 'stddev_epi', '-dpng')
montage1 = createMontage(tSNR_3D, 5, 1, 'tSNR (whole image)', 'hot');
print(montage1.f, 'tsnr_whole', '-dpng')
montage4 = createMontage(tSNR_3D_masked, 5, 1, 'tSNR (brain)', 'hot');
print(montage4.f, 'tsnr_brain', '-dpng')
figmask = displayMaskContour(F3D_mean, mask_3D, 0, 3);
print(figmask, 'mask_contour', '-dpng')



% F2D_mean_masked = zeros(Ni*Nj*Nk, 1);
% F2D_mean_masked(I_mask, :) = F2D_mean(I_mask, :);
% F3D_mean_masked = reshape(F2D_mean_masked, Ni, Nj, Nk);
% montage5 = createMontage(F3D_mean_masked, 5, 1, 'Mean EPI (brain)');

% F2D_var_masked = zeros(Ni*Nj*Nk, 1);
% F2D_var_masked(I_mask, :) = F2D_var(I_mask, :);
% F3D_var_masked = reshape(F2D_var_masked, Ni, Nj, Nk);
% montage6 = createMontage(F3D_var_masked, 5, 1, 'Variance (brain)');

% F2D_stddev_masked = zeros(Ni*Nj*Nk, 1);
% F2D_stddev_masked(I_mask, :) = F2D_stddev(I_mask, :);
% F3D_stddev_masked = reshape(F2D_stddev_masked, Ni, Nj, Nk);
% montage7 = createMontage(F3D_stddev_masked, 5, 1, 'Standard deviation (brain)');

% Create html log file
log_nr = 5;
log_name = [subject '_log' num2str(log_nr) '.html'];
fid = fopen(log_name,'a');
fprintf(fid, '<H2>Log</H2>');
fprintf(fid, ['\n<BR>Subject:  ' subject]);
t = datestr(datetime('now'));
fprintf(fid, ['\n<BR>Date/time:  ' t]);

fprintf(fid, '<H2>Imaging info</H2>');
fprintf(fid, ['\nVolumes:  ' num2str(Nt)]);
fprintf(fid, ['\n<BR>Voxels (x,y,z):  ' num2str(Ni) ', ' num2str(Nj) ', ' num2str(Nk)]);

fprintf(fid, '<H2>Timeseries summary</H2>');
fprintf(fid, '\n<TABLE><TR><TD><img src="timeseries_summary.png" alt="no picture" width=700 height=600></TD></TR></TABLE>' );

fprintf(fid, '<H2>QC Metrics</H2>');
fprintf(fid, ['\nFD threshold (mm):  ' num2str(FD_threshold)]);
fprintf(fid, ['\n<BR>FD outliers:  ' num2str(numel(FD_measures.FD_outliers_ind))]);
fprintf(fid, ['\n<BR>Total FD:  ' num2str(FD_measures.FD_sum)]);
fprintf(fid, ['\n<BR>Mean FD:  ' num2str(FD_measures.FD_mean)]);
fprintf(fid, ['\n<BR>Mean Zscore:  ' num2str(Zstat_mean)]);
fprintf(fid, ['\n<BR>tSNR (brain):  ' num2str(tSNR_brain)]);
fprintf(fid, ['\n<BR>tSNR (GM):  ' num2str(tSNR_GM)]);
fprintf(fid, ['\n<BR>tSNR (WM):  ' num2str(tSNR_WM)]);
fprintf(fid, ['\n<BR>tSNR (CSF):  ' num2str(tSNR_CSF)]);

fprintf(fid, '<H2>QC brain images</H2>');
fprintf(fid, '\n<TABLE><TR><TD><img src="mean_epi.png" alt="no picture" width=700 height=600></TD></TR></TABLE>' );
fprintf(fid, '\n<TABLE><TR><TD><img src="stddev_epi.png" alt="no picture" width=700 height=600></TD></TR></TABLE>' );
fprintf(fid, '\n<TABLE><TR><TD><img src="tsnr_whole.png" alt="no picture" width=700 height=600></TD></TR></TABLE>' );
fprintf(fid, '\n<TABLE><TR><TD><img src="tsnr_brain.png" alt="no picture" width=700 height=600></TD></TR></TABLE>' );
fprintf(fid, '\n<TABLE><TR><TD><img src="mask_contour.png" alt="no picture" width=700 height=700></TD></TR></TABLE>' );
fclose(fid);


























