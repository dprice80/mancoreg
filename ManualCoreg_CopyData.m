clear all
close all
clc

addpath /imaging/camcan/sandbox/projects/QueryFunction/QueryFun_v1/

DAT.SessionList = {
    'MRI' '/imaging/camcan/cc700/mri/pipeline/release003/data/aamod_coreg_t2_structural_00001/<MRIID>/structurals/rs*.nii'
    };
DAT = CCQuery_CheckFiles(DAT);

%% Copy Data
DAT.CopyData.DestRoot = '/imaging/camcan/sandbox/projects/megcoregister/data/anatomical_coreg/release003/';
DAT.CopyData.FileList = find(all([DAT.FileCheck],2))'; % find the indexes of files that exist for all file checks
DAT.CopyData.DestSubDirs = {
    'MRI' 'anatomicals'
    }; 
DAT.CopyData.SoftLink = 1;
DAT.CopyData.Overwrite = 0;
DAT = CCQuery_CopyData(DAT);