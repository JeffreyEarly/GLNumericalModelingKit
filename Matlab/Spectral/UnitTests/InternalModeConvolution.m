[rhoFunc, N2Func, zIn] = InternalModes.StratificationProfileWithName('exponential');
N = 32;
k = 0;

z = linspace(min(zIn),max(zIn),1024).';
im = InternalModesSpectral(rhoFunc,zIn,z,33);

z = im.GaussQuadraturePointsForModesAtWavenumber(N,0);
im = InternalModes(rhoFunc,zIn,z,33,'nModes',N);


%%%%%
% NOTE: Two G modes do not project onto a G mode, they project onto F
% modes. Think odd-even sine-cosine,

[F,G] = im.ModesAtWavenumber(0);

% add the barotropic mode
F = circshift(F,1,2);
F(:,1) = 1;

n_eff = N-2;
M = G(:,1:n_eff);

n_eff = N;
M = F;

fn = zeros(n_eff,1);
fn(2) = 2;
f = M*fn;

gn = zeros(n_eff,1);
gn(2) = 2;
g = M*gn;

a_out = M\(f.*g);

figure
plot(f.*g,z)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% We are going to compute the interaction sum between F, F and F modes.
% The ordering will be such that the third index is the destination.

C = ConstructInteractionMatrix(F,F,F);

fn = randn(N,1);
gn = randn(N,1);

f = F*fn;
g = F*gn;

hn_spatial = F\(f.*g);

hn_spectral = zeros(N,1);
for i=1:N
   hn_spectral(i) = (fn.')*C(:,:,i)*gn ;
end

%%%%%%%%%%%%%%%%
% cool plot
figure, pcolor(log10(abs(C(2:end,2:end,5))))
colorbar('eastoutside')