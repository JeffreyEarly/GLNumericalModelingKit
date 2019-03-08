knot_type = 6;

% The data
t = (0:10)';
x = [0;0;0;0;0;0;2;4;6;8;10];

% Problem constants
tq = linspace(0,max(t),1000)';
K = 3;
S = K-1;

group = struct('left',[],'right',[],'K',[]);

if knot_type == 1 % 1. Natural knots
    t_knot = NaturalKnotsForSpline(t,K,1);
    constraints = struct('t',[],'D',[]);
elseif knot_type == 2
    % This works well. It shows three distinct accelerations, which is not
    % how we grouped things, but it's assuming they're outside our
    % observation window.
    t_knot = [0;4;4.5;5.5;6.;10]; %%% Might be my favorite!!??
    constraints = struct('t',[2.5;2.5;7.5],'D',[1;2;2]);
elseif knot_type == 3
    % We don't like this because it supposes that a group of points that we
    % said wasn't accelerating, was in fact, accelerating, it's just that
    % we missed it.
    t_knot = [0;3.5;4.5;5.5;6.5;10];
    constraints = struct('t',[2.5;2.5;7.5],'D',[1;2;2]);
elseif knot_type == 4
    % 
    t_knot = [0;5;5;10];
    constraints = struct('t',[2.5;2.5;7.5],'D',[1;2;2]);
elseif knot_type == 5
    % The captures the essense of what we want just fine, but it has a
    % discontinuous velocity.
    t_knot = [0;4;4;6;6;10];
       constraints = struct('t',[2.5;2.5;7.5],'D',[1;2;2]);
elseif knot_type == 6
    % The captures the essense of what we want just fine, but it has a
    % discontinuous velocity.
    t_knot = [0;4.5;4.5;5.5;5.5;10];
    constraints = struct('t',[2.5;2.5;7.5],'D',[1;2;2]);
elseif knot_type == 6
    % The captures the essense of what we want just fine, but it has a
    % discontinuous velocity.
    t_knot = [0;4.5;4.5;5.5;5.5;10];
    constraints = struct('t',[2.5;2.5;7.5],'D',[1;2;2]);
end


% The 'natural' knots are the ideal situation, where you assume the
% underlying data is C^\infnty, and so you interpolate to the best of your
% ability. With that as the context, the transition knots are the obvious
% choice. The argument against them, is just that you don't think that the
% underlying data is necessary C^3 everywhere. In that case, you want
% discontinuities--probably at the midway knots.

f = ConstrainedSpline(t,x,K,t_knot,NormalDistribution(1),constraints);

fig1 = figure;

subplot(4,1,1)
plot(tq,f(tq),'LineWidth',2), hold on
scatter(t,x,6^2,'filled')
xlim([t(1) t(end)])
ylim([-1 11])
% set(gca, 'XTick', []);
set(gca, 'XTickLabel', []);
set(gca, 'YTick', []);
vlines(t_knot,'g--')
ylabel('position');

subplot(4,1,2)
plot(tq,f(tq,1),'LineWidth',2), hold on
xlim([t(1) t(end)])
% set(gca, 'XTick', []);
set(gca, 'XTickLabel', []);
set(gca, 'YTick', []);
ylabel(sprintf('$1^{st}$-deriv'),'Interpreter','LaTex');
vlines(t_knot,'g--')

if K > 2
    subplot(4,1,3)
    plot(tq,f(tq,2),'LineWidth',2), hold on
    xlim([t(1) t(end)])
    % set(gca, 'XTick', []);
    set(gca, 'XTickLabel', []);
    set(gca, 'YTick', []);
    ylabel(sprintf('$2^{nd}$-deriv'),'Interpreter','LaTex');
    vlines(t_knot,'g--')
end

% subplot(4,1,4)
% plot(tq,Bq(:,:,1),'LineWidth',2), hold on
% xlim([t(1) t(end)])
% set(gca, 'YTick', []);
% ylabel('b-splines', 'FontSize', figure_axis_label_size, 'FontName', figure_font);
% vlines(t_knot,'g--')

if length(group.K) > 0
    showGroupsOnFigure(group, fig1, t, Helvetica, 9)
end