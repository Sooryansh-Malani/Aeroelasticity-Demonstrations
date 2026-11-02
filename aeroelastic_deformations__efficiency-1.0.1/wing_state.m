function [s, dy, thx, thz, v_div] = wing_state(parm, N, v, sweep_deg, ea_frac)
%WING_STATE  Static aeroelastic state of the wing at one speed and one sweep angle.
%
%   [s,dy,thx,thz,v_div] = wing_state(parm,N,v,sweep_deg,ea_frac)
%
%   Inputs
%     parm       parameter struct from the main script
%     N          number of strips/elements (e.g. 100)
%     v          airspeed [m/s]  (must be below the divergence speed)
%     sweep_deg  sweep angle [deg]
%     ea_frac    elastic axis position as a fraction of chord from the LE
%                (optional, default 0.3)
%
%   Outputs (row vectors, one entry per node, node N = wing tip)
%     s      position of each node along the elastic axis [m]
%     dy     vertical deflection [m]
%     thx    twist about the elastic axis [rad]
%     thz    bending slope [rad]
%     v_div  divergence speed [m/s]
%
%   For a value at any other point use interpolation, e.g.
%     interp1([0 s],[0 dy],s_query)
if nargin < 5 || isempty(ea_frac), ea_frac = 0.3; end
if ~isfield(parm,'alpha'), parm.alpha = deg2rad(5); end

parm.swept         = deg2rad(sweep_deg);
parm.v             = v;
parm.elastic_chord = ea_frac*parm.chord_length;
parm.e             = abs(parm.elastic_chord - parm.ac_chord);

[dy, thz, thx, ~, v_div] = deformation(N, parm);

s   = (1:N)*parm.span/N;
dy  = dy(:).';
thz = thz(:).';
thx = thx(:).';
end