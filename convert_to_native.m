function fid = convert_to_native(coreg_im, native_im, fidlocs)
% Convert coregistered fiducials to native space. 
% coreg_im = coregistered non resliced image
% native_im = native space image (raw nii file) 
% fidlocs = the output from ManualCoregFun

% MNI space T
hdr = load_nii_hdr(coreg_im); % coregistered non resliced image
T = [hdr.hist.srow_x; hdr.hist.srow_y; hdr.hist.srow_z; [0 0 0 1]];
% Native T
hdr = load_nii_hdr(native_im); % native space image
Tnative = [hdr.hist.srow_x; hdr.hist.srow_y; hdr.hist.srow_z; [0 0 0 1]];

F = load([fileparts(fidlocs), '/FiducialLocs.mat']);
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
    
end