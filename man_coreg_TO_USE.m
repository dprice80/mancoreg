%% Manual registration of fiducials for M/EEG coregistration
%%Uses query function
%%run in matlab 2015, errors likely later versions
%%Thanks to Darren Price for creation Sept 2018

addpath /imaging/dp01/toolboxes/query_generic/
addpath /imaging/camcan/cc700/mri/pipeline/release003/megcoregister/Scripts/nifti_analyse
addpath /imaging/camcan/cc700/mri/pipeline/release003/megcoregister/Scripts/

q = queryfs('/imaging/dn01/Documents/All_final_analysis_July_2018_onwards/HMM_Roni_collab_v1.0/idlist-14-9-18.csv');
q.addsearchpath('T1','/imaging/camcan/ccfrail/mri/pipeline/release001/func/aamod_coreg_extended_1_00001/<ID>/structurals/sMR*.nii');
q.selectfirstfile = true;
q.checkfiles

for ii = q.allexist
    
    mkdir(sprintf('/imaging/dn01/Documents/All_final_analysis_July_2018_onwards/HMM_Roni_collab_v1.0/manualcoreg/%s'), q.ID.ID{ii})
    unix(sprintf('cp %s /imaging/dn01/Documents/All_final_analysis_July_2018_onwards/HMM_Roni_collab_v1.0/manualcoreg/%s', q.filenames.T1{ii}, q.ID.ID{ii}))
    
    flags.which = 1;
    [fl, fil, ext] = fileparts(q.filenames.T1{ii});
    spm_reslice({'/imaging/dn01/Documents/All_final_analysis_July_2018_onwards/HMM_Roni_collab_v1.0/c1MNI152_T1_1mm_brain.nii', ['/imaging/dn01/Documents/All_final_analysis_July_2018_onwards/HMM_Roni_collab_v1.0/manualcoreg/' q.ID.ID{ii} '/' fil ext]}, flags)
    
end



%%

q2 = [];
q2 = queryfs('/imaging/dn01/Documents/All_final_analysis_July_2018_onwards/HMM_Roni_collab_v1.0/idlist-14-9-18.csv');
q2.addsearchpath('T1rs','/imaging/dn01/Documents/All_final_analysis_July_2018_onwards/HMM_Roni_collab_v1.0/manualcoreg/<ID>/rs*');
q2.addsearchpath('complete','/imaging/dn01/Documents/All_final_analysis_July_2018_onwards/HMM_Roni_collab_v1.0/manualcoreg/<ID>/Fiducial*');
q2.addsearchpath('T1','/imaging/camcan/ccfrail/mri/pipeline/release001/func/aamod_coreg_extended_1_00001/<ID>/structurals/sMR*.nii');
q2.addsearchpath('T1native','/imaging/camcan/ccfrail/mri/pipeline/release001/data/aamod_convert_structural_00001/<ID>/structurals/sMR*.nii');
q2.checkfiles

T1rsind = q2.fileindex.T1rs;
T1rsind(68) = false;

for ii = find(T1rsind & ~q2.fileindex.complete)'
    ManualCoregFun(q2.filenames.T1rs{ii});
end



%% This will save coordinates in voxel and mm units
for ii = find(T1rsind)'
    % MNI space T
    hdr = load_nii_hdr(q2.filenames.T1{ii});
    T = [hdr.hist.srow_x; hdr.hist.srow_y; hdr.hist.srow_z; [0 0 0 1]];
    % Native T
    hdr = load_nii_hdr(q2.filenames.T1native{ii});
    Tnative = [hdr.hist.srow_x; hdr.hist.srow_y; hdr.hist.srow_z; [0 0 0 1]];
    
    F = load([fileparts(q2.filenames.T1rs{ii}), '/FiducialLocs.mat']);
    if isfield(F, 'Mmm')
        natvox = [F.Mmm [1 1 1]'] * inv(T)';
        
        fid = [];
        fid.native.vox.lpa = natvox(1,1:3); % native space voxel coords
        fid.native.vox.rpa = natvox(2,1:3);
        fid.native.vox.nas = natvox(3,1:3);
        fid.mni.mm.lpa     = F.Mmm(1,1:3); % save mni space mm coords
        fid.mni.mm.rpa     = F.Mmm(2,1:3);
        fid.mni.mm.nas     = F.Mmm(3,1:3);
        natmm = natvox * Tnative'; % transform voxel coords to native space mm
        fid.native.mm.lpa = natmm(1,1:3);
        fid.native.mm.rpa = natmm(2,1:3);
        fid.native.mm.nas = natmm(3,1:3);
        save([fileparts(q2.filenames.T1rs{ii}) '/FiducialLocs.mat'], 'fid')
    end
end

%% 


q2 = [];
q2 = queryfs('/imaging/dn01/Documents/All_final_analysis_July_2018_onwards/HMM_Roni_collab_v1.0/idlist-14-9-18.csv');
q2.addsearchpath('T1rs','/imaging/dn01/Documents/All_final_analysis_July_2018_onwards/HMM_Roni_collab_v1.0/manualcoreg/<ID>/rs*');
q2.addsearchpath('complete','/imaging/dn01/Documents/All_final_analysis_July_2018_onwards/HMM_Roni_collab_v1.0/manualcoreg/<ID>/Fiducial*');
q2.addsearchpath('T1','/imaging/camcan/ccfrail/mri/pipeline/release001/func/aamod_coreg_extended_1_00001/<ID>/structurals/sMR*.nii');
q2.addsearchpath('T1native','/imaging/camcan/ccfrail/mri/pipeline/release001/data/aamod_convert_structural_00001/<ID>/structurals/sMR*.nii');
q2.checkfiles


colormap gray
for ii = find(T1rsind)'
    nii = load_untouch_nii(q2.filenames.T1{ii});
    load([fileparts(q2.filenames.T1rs{ii}) '/FiducialLocs.mat']);
    figure(1); clf
    imagesc(squeeze(nii.img(:,round(fid.native.vox.lpa(2)),:)))
    set(gca,'YDir','normal')
    hold on
    plot(fid.native.vox.lpa(3), fid.native.vox.lpa(1), 'Marker', '+', 'MarkerSize', 30)
    pause(1)
end

%%
% Load resliced coregistered image
view_nii(load_nii(q2.filenames.T1rs{109}));

% view slice from non-resliced
nii = load_untouch_nii(q2.filenames.T1{109});

figure
imagesc(nii.img(:,:,120))

