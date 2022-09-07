
%% create b0 scan files (TOPPE)

addpath ~/github/HarmonizedMRI/Calibration/b0/GE   % b04ge.m

sys = toppe.systemspecs('maxSlew', 15, 'slewUnit', 'Gauss/cm/ms', ...
    'maxGrad', 5, 'gradUnit', 'Gauss/cm', ...
    'myrfdel', 58, ... % (us) best if multiple of 4us
    'daqdel', 60, ...  % (us) best if multiple of 4us
    'gradient', 'xrm', ... % xrm: MR750; hrmb: UHP
    'timessi', 200);    % us


% isotropic FOV and matrix
FOV = 24 * [1 1 1];  % cm  
N = 120 * [1 1 1]; 

deltaTE = [0 2.3];          % (ms) TE shift for B0 mapping
flip = 5;                   % excitation angle (degrees)

b04ge(sys, N, FOV, flip, deltaTE, ...
    'tbw', 12, ...
    'rfDur', 0.75, ...   % ms
    'ftype', 'min', ...
    'slabThick', FOV(3) * 0.9, ...
    'rfSpoilSeed', 117, ...
    'nCyclesSpoil', 1, ... % spoiling gradient cycles across one voxel
    'fatsat', false);

system('mv b0.tar localizer.tar');
