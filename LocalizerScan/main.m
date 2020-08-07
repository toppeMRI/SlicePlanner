function main
% create 3D spin-warp toppev3 sequence
% isotropic matrix and voxel dimensions
% acquires B0 field map as well

% addpath ~/github/toppe/

%% set system limit structs
% NB! If passing 'sys' to writemod.m, 'maxGrad' MUST match the physical system limit -- since gradients are scaled relative to this.
% 'maxSlew' can always be a design choice, i.e., it can be at or below the physical system limit.
% Here we will therefore use 'sys' when designing waveforms, and 'sysGE' when writing them to .mod files with writemod.m.
mxs = 10.0;    % max slew (G/cm/msec). Go easy to minimize PNS.
seq.sys = toppe.systemspecs('maxSlew', mxs, 'slewUnit', 'Gauss/cm/ms', 'maxGrad', 5, 'gradUnit', 'Gauss/cm');

%% Create modules.txt
% Entries are tab-separated.
modFileText = ['' ...
'Total number of unique cores\n' ...
'2\n' ...
'fname	duration(us)	hasRF?	hasDAQ?\n' ...
'readout.mod	0	0	1\n' ...
'tipdown.mod	0	1	0' ];
fid = fopen('modules.txt', 'wt');
fprintf(fid, modFileText);
fclose(fid);

%% Sequence parameters
seq.n = 120;                     % matrix size (isotropic)
seq.fov = 24;                    % fov (cm) (isotropic)
seq.oprbw = 125/4;               % kHz
dx = seq.fov/seq.n;              % voxel size (cm)
seq.deltaTE = 2.3;               % (ms) TE shift for B0 mapping
seq.nCycleSpoil = 1;            % readout spoiler gradient

seq.rf.flip = 5;                   % excitation angle (degrees)
seq.rf.slabThick = seq.fov*0.9;    % cm
seq.rf.tbw = 12;                   % time-bandwidth product of SLR pulse 
seq.rf.dur = 0.75;                 % RF pulse duration (msec)
seq.rf.ftype = 'min';              % a good option for 3D imaging is 'min'; otherwise 'ls' is common
seq.rf.nCycleSpoil = seq.n*seq.nCycleSpoil;

%% Create sequence modules
[rf,gex] = toppe.utils.rf.makeslr(seq.rf.flip, seq.rf.slabThick, seq.rf.tbw, seq.rf.dur, seq.rf.nCycleSpoil, ...
	'ftype',  seq.rf.ftype, ...
	'spoilDerate', 0.5, ...
	'system', seq.sys, ...
	'ofname', 'tipdown.mod');

% Remember: makegre creates y and z phase-encodes based on isotropic resolution
[gx,gy,gz] = toppe.utils.makegre(seq.fov, seq.n, seq.fov/seq.n, ...
	'system',  seq.sys, ...
	'ncycles', seq.nCycleSpoil, ...
	'ofname', 'readout.mod', ...
	'slewDerate', 0.5, ...
	'oprbw',   seq.oprbw);

%% Create scanloop.txt
rfphs = 0;              % radians
rf_spoil_seed_cnt = 0;
seq.rf_spoil_seed = 117;

ny = seq.n;
nz = seq.n;

fprintf('Writing scanloop.txt \n');
toppe.write2loop('setup', 'version', 3);
for iz = -0:nz           % iz < 1 are discarded acquisitions (to reach steady state)
	if ~mod(iz,10)
		fprintf([repmat('\b',1,30) sprintf('\t%d of %d', iz, nz)]);
	end
	if iz > 0
		dabmode = 'on';
		peOn = true;
	else
		dabmode = 'off';
		peOn = false;
	end

	for iy = 1:ny
		% y/z phase-encode amplitudes, scaled to (-1,1)
		a_gy = peOn * ((iy-1+0.5)-ny/2)/(ny/2);
		a_gz = peOn * ((iz-1+0.5)-nz/2)/(nz/2);

		for iim = 1:2
			% TE shift for B0 field mapping
			if iim == 1
				textra1 = 0;
				textra2 = seq.deltaTE; 
			else
				textra1 = seq.deltaTE; 
				textra2 = 0;
			end

			% rf excitation
	  		toppe.write2loop('tipdown.mod', ...
				'RFamplitude', 1.0, ...
				'RFphase', rfphs, ...
				'textra', textra1);

		 	% readout. Data is stored in 'slice', 'echo', and 'view' indeces.
			toppe.write2loop('readout.mod', ...
				'Gamplitude', [1 a_gy a_gz]', ...
				'DAQphase', rfphs, ...
				'textra', textra2, ...
				'slice', max(iz,1), 'echo', iim, 'view', iy, ...
				'dabmode', dabmode);

			% update rf phase (RF spoiling)
			rfphs = rfphs + (seq.rf_spoil_seed/180 * pi)*rf_spoil_seed_cnt ;  % radians
			rf_spoil_seed_cnt = rf_spoil_seed_cnt + 1;
		end
	end
end
fprintf('\n');
toppe.write2loop('finish');

%% Save parameters and create tar file.
save seq seq
system('tar czfP toppev3,localizer.tgz main.m seq.mat *.mod modules.txt scanloop.txt');
system('rm seq.mat');

fprintf('Scan time for 3D sequence: %.2f min\n', toppe.getscantime/60);
fprintf('Copy toppev3,localizer.tgz to /usr/g/bin/ on scanner, untar, and scan with toppev3.\n');

return;
