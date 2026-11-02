function res = plot_wing_deformations(parm, v, sweeps, ea_frac)
%PLOT_WING_DEFORMATIONS  3-D picture of the deformed wing at several sweep angles.
%
%   res = plot_wing_deformations(parm, v, sweeps, ea_frac)
%
%   parm     parameter struct from the main script
%   v        airspeed [m/s]               (default 30)
%   sweeps   sweep angles [deg]           (default [-10 0 20])
%   ea_frac  elastic axis, fraction of chord from the LE (default 0.3)
%
%   Every wing is solved with N = 100 strips, but only every 10th section is
%   drawn (plus the clamped root). Each drawn section is two points, leading
%   and trailing edge. They move by
%       dz = delta_y + (distance ahead of the elastic axis) * d_alpha
%   so the LE rises and the TE drops for nose-up twist.
%
%   res(k) holds sweep, s, dy, thx, thz, v_div and the tip values.
%
%   Axes: x = streamwise (flow toward +x), y = span, z = vertical.

if nargin < 2 || isempty(v),       v = 30;               end
if nargin < 3 || isempty(sweeps),  sweeps = [-10 0 20];  end
if nargin < 4 || isempty(ea_frac), ea_frac = 0.3;        end
if ~isfield(parm,'alpha'), parm.alpha = deg2rad(5); end

%% Settings
N           = 100;      % strips used for the solution
step        = 10;       % draw every 10th section
scale_dy    = 1;        % deflection exaggeration (1 = true scale)
scale_twist = 10;       % twist exaggeration (the real twist is tiny, see below)
use_bending_slope = true;   % true : d_alpha = cos(L)*thx - sin(L)*thz  (as in the model)
                            % false: d_alpha = thx only

c      = parm.chord_length;
idx    = [0, step:step:N];          % section 0 = clamped root
nS     = numel(sweeps);
colors = lines(nS);

X = cell(1,nS);  Y = cell(1,nS);  Z = cell(1,nS);  Xea = cell(1,nS);
Yea = cell(1,nS); Zea = cell(1,nS);
res = struct('sweep',{},'s',{},'dy',{},'thx',{},'thz',{},'v_div',{}, ...
             'tip_dy',{},'tip_twist_deg',{});

%% Solve and build the geometry
for k = 1:nS
    sw  = sweeps(k);
    lam = deg2rad(sw);

    [s,dy,thx,thz,v_div] = wing_state(parm,N,v,sw,ea_frac);
    if v >= v_div
        warning('Sweep %g deg: V = %g m/s is above the divergence speed (%.1f m/s). Result is not meaningful.', ...
                sw, v, v_div);
    end

    if use_bending_slope
        dalpha = cos(lam)*thx - sin(lam)*thz;
    else
        dalpha = thx;
    end

    s_all  = [0 s];
    dy_all = [0 dy];
    da_all = [0 dalpha];
    sel    = idx + 1;

    sp   = s_all(sel);
    Yea_ = sp*cos(lam);                 % span position
    Xea_ = sp*sin(lam);                 % streamwise position of the elastic axis
    Xle  = Xea_ - ea_frac*c;
    Xte  = Xea_ + (1-ea_frac)*c;

    zea  = scale_dy*dy_all(sel);
    Zle  = zea + ea_frac*c*scale_twist*da_all(sel);
    Zte  = zea - (1-ea_frac)*c*scale_twist*da_all(sel);

    X{k}   = [Xle(:) Xte(:)];
    Y{k}   = [Yea_(:) Yea_(:)];
    Z{k}   = [Zle(:) Zte(:)];
    Xea{k} = Xea_;  Yea{k} = Yea_;  Zea{k} = zea;

    res(k).sweep         = sw;
    res(k).s             = s;
    res(k).dy            = dy;
    res(k).thx           = thx;
    res(k).thz           = thz;
    res(k).v_div         = v_div;
    res(k).tip_dy        = dy(end);
    res(k).tip_twist_deg = rad2deg(thx(end));
end

%% Common axis limits
allX = cellfun(@(a) [min(a(:)) max(a(:))], X, 'UniformOutput', false);  allX = vertcat(allX{:});
allY = cellfun(@(a) [min(a(:)) max(a(:))], Y, 'UniformOutput', false);  allY = vertcat(allY{:});
allZ = cellfun(@(a) [min(a(:)) max(a(:))], Z, 'UniformOutput', false);  allZ = vertcat(allZ{:});
xl = [min(allX(:,1)) max(allX(:,2))];
yl = [0 max(allY(:,2))];
zl = [min([allZ(:);0]) max([allZ(:);0])];
pad = 0.05*max([diff(xl) diff(yl) diff(zl)]);
xl = xl + [-pad pad];  zl = zl + [-pad pad];

%% Plot
figure('Color','w','Name','Deformed wing');
tiledlayout(1,nS,'TileSpacing','compact','Padding','compact');
for k = 1:nS
    ax = nexttile; hold(ax,'on');

    h0 = surf(ax,X{k},Y{k},zeros(size(Z{k})),'FaceColor',[0.8 0.8 0.8], ...
              'FaceAlpha',0.35,'EdgeColor',[0.55 0.55 0.55]);
    h1 = surf(ax,X{k},Y{k},Z{k},'FaceColor',colors(k,:),'FaceAlpha',0.75, ...
              'EdgeColor','k');
    h2 = plot3(ax,Xea{k},Yea{k},Zea{k},'r-','LineWidth',1.5);
    plot3(ax,X{k}(:),Y{k}(:),Z{k}(:),'ko','MarkerFaceColor','k', ...
          'MarkerSize',3,'HandleVisibility','off');

    if res(k).v_div > 1000
        vds = '> 1000 m/s';
    else
        vds = sprintf('%.0f m/s',res(k).v_div);
    end
    title(ax,{sprintf('\\Lambda = %g^\\circ,  V_D = %s',res(k).sweep,vds), ...
              sprintf('tip: \\delta_y = %.1f mm,  \\theta_x = %.2f^\\circ', ...
                      1e3*res(k).tip_dy,res(k).tip_twist_deg)}, ...
          'Color','k','FontWeight','bold');
    xlabel(ax,'x - streamwise [m]','Color','k','FontWeight','bold');
    ylabel(ax,'y - span [m]','Color','k','FontWeight','bold');
    zlabel(ax,'z - vertical [m]','Color','k','FontWeight','bold');
    set(ax,'FontWeight','bold','XColor','k','YColor','k','ZColor','k');
    xlim(ax,xl); ylim(ax,yl); zlim(ax,zl);
    daspect(ax,[1 1 1]);
    view(ax,-37.5,25); grid(ax,'on'); box(ax,'on');

    if k == 1
        legend(ax,[h0 h1 h2],{'undeformed','deformed','elastic axis'},'Location','northwest', ...
               'TextColor',[1 1 0],'FontWeight','bold');
    end
end
sgtitle(sprintf('Deformed wing at V = %g m/s, \\alpha = %g^\\circ   (deflection x%g, twist x%g)', ...
        v, rad2deg(parm.alpha), scale_dy, scale_twist), 'Color','k','FontWeight','bold');
end