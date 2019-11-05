function main

% create 3D spin-warp toppev3 sequence
% isotropic matrix and voxel dimensions

% addpath ~/github/toppe/

%% set system limit structs
% NB! If passing 'sys' to writemod.m, 'maxGrad' MUST match the physical system limit -- since gradients are scaled relative to this.
% 'maxSlew' can always be a design choice, i.e., it can be at or below the physical system limit.
% Here we will therefore use 'sys' when designing waveforms, and 'sysGE' when writing them to .mod files with writemod.m.
mxs = 12.0;    % max slew [G/cm/msec]. Go easy to minimize PNS.
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
seq.fov = 24;                    % in-plane fov (cm)
seq.oprbw = 125/4;               % kHz
dx = seq.fov/seq.n;              % voxel size (cm)

seq.ncyclesSpoil = 1;              % readout spoiler gradient

seq.rf.flip = 5;                   % excitation angle (degrees)
seq.rf.slabThick = seq.fov*0.9;    % cm
seq.rf.tbw = 12;                    % time-bandwidth product of SLR pulse 
seq.rf.dur = 0.75;                    % RF pulse duration (msec)
seq.rf.ftype = 'min';              % a good option for 3D imaging is 'min'; otherwise 'ls' is common
seq.rf.ncyclesSpoil = seq.n*seq.ncyclesSpoil;

%% Create sequence modules
[rf,gex] = toppe.utils.rf.makeslr(seq.rf.flip, seq.rf.slabThick, seq.rf.tbw, seq.rf.dur, seq.rf.ncyclesSpoil, ...
                       'ftype', seq.rf.ftype, 'system', seq.sys, 'ofname', 'tipdown.mod');

% Remember: makegre creates y and z phase-encodes based on isotropic resolution
[gx,gy,gz] = toppe.utils.makegre(seq.fov, seq.n, seq.fov/seq.n, 'system', seq.sys, 'ncycles', seq.ncyclesSpoil, 'oprbw', seq.oprbw);

%% Create scanloop.txt
rfphs = 0;              % radians
rf_spoil_seed_cnt = 0;
rf_spoil_seed = 117;

ny = seq.n;
nz = seq.n;

fprintf('Writing scanloop.txt for B0 sequence\n');
toppe.write2loop('setup', 'version', 3);
for iz = -0:nz           % iz < 1 are discarded acquisitions (to reach steady state)
	if ~mod(iz,10)
		fprintf([repmat('\b',1,20) sprintf('%d of %d', iz, nz)]);
	end

	for iy = 1:ny
		for iim = 1:1
			switch iim
				case 1
					textra1 = 0;
					textra2 = 0; %2.3 + tdelay;     % msec
				case 2
					textra1 = 2.3;     % msec
					textra2 = 0 + tdelay;
			end

         if iz > 0
            a_gy = ((iy-1+0.5)-ny/2)/(ny/2);    % y phase-encode amplitude, scaled to (-1,1) range
        		a_gz = ((iz-1+0.5)-nz/2)/(nz/2);    % z phase-encode amplitude, scaled to (-1,1) range
				dabmode = 'on';
         else
            a_gy = 0;
         	a_gz = 0; 
				dabmode = 'off';
         end

	   	% rf excitation
	  		toppe.write2loop('tipdown.mod', 'RFphase', rfphs, 'RFamplitude', 1.0, 'textra', textra1);

		 	% readout. Data is stored in 'slice', 'echo', and 'view' indeces.
   		toppe.write2loop('readout.mod', 'DAQphase', rfphs, ...
				'slice', max(iz,1), 'echo', iim, 'view', iy, ...
				'dabmode', dabmode, 'Gamplitude', [1 a_gy a_gz]', 'textra', textra2);

		   % update rf phase (RF spoiling)
			rfphs = rfphs + (rf_spoil_seed/180 * pi)*rf_spoil_seed_cnt ;  % radians
			rf_spoil_seed_cnt = rf_spoil_seed_cnt + 1;
      end
	end
end
fprintf('\n');
toppe.write2loop('finish');

%% Save parameters and create tar file. We can reuse modules.txt from above.
system('mkdir -p tar');
cd tar
save seq seq
system('tar cf toppev3,localizer.tar ../main.m seq.mat ../*.mod ../modules.txt ../scanloop.txt');
cd ..

fprintf('Scan time for 3D sequence: %.2f min\n', toppe.getscantime/60);

return;

