% Example to run manual coreg on some camcan data
% Must use matlab 2015a (doesn't work in 2018a)

clear all

restoredefaultpath

addpath /imaging/dp01/toolboxes/mancoreg/nifti_analyse
addpath /imaging/dp01/toolboxes/mancoreg
mkdir /imaging/dp01/toolboxes/mancoreg_example

!rsync -a /imaging/camcan/cc700/mri/pipeline/release004/data/aamod_convert_structural_00001/CC110033/structurals/sMR10033_CC110033-0003-00001-000192-01.nii /imaging/dp01/toolboxes/mancoreg_example

example_file = '/imaging/dp01/toolboxes/mancoreg_example/sMR10033_CC110033-0003-00001-000192-01.nii';


%% 

ManualCoregFun(example_file);

load('sMR10033_CC110033-0003-00001-000192-01_fids.mat')


