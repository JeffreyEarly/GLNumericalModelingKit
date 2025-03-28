N = 1024;
Nz = 500;
x = rand(N,N,Nz);
nLoops = 4;

nCores = (1:16).';
wallTime = zeros(size(nCores));

for iCore=1:length(nCores)
    dft = RealToComplexTransform(size(x),dims=[1 2],nCores=nCores(iCore),planner="estimate");
    tic
    for i=1:nLoops
        y2 = dft.transformForward(x);
    end
    wallTime(iCore) = toc;
end

figure, scatter(nCores,wallTime,4^2,'filled')