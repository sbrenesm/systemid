[b,a]=cheby1(3,0.5,[2*0.15 2*0.3]);
Nest = 1000;

    u0e = randn(Nest,1);
    nye = normrnd(0,0.5,Nest,1);
    y0e = filter(b,a,u0e);
    ye = y0e+nye;

error = experiment1(a,b);


function [error]= experiment1(a,b)
    %size of I
    gsize = 100;
    %estimation dataset
    Nest = 1000;
    u0e = randn(Nest,1);
    nye = normrnd(0,0.5,Nest,1);
    y0e = filter(b,a,u0e);
    ye = y0e+nye;
    
    %evaluation dataset
    Nval = 10000;
    u0v = randn(Nval,1);
    nyv = normrnd(0,0.5,Nval,1);
    y0v = filter(b,a,u0v);
    yv = y0v+nyv;
    
    %create matrices
    Umate = zeros(Nest,gsize);
    Umatv = zeros(Nval,gsize);
    g =zeros(gsize,gsize);
    error = zeros(gsize,1);
    errorak = zeros(gsize,1);
    %fill matrices
    for i =1:gsize
        Umate(i:Nest,i) = u0e(1:Nest-i+1);
        Umatv(i:Nval,i) = u0v(1:Nval-i+1);
    end
    %solve for the evaluation dataset and compute error
    for i =1:100
        g(1:i,i) = Umate(:,1:i)\ye;  
        error(i) = (transpose(yv-mtimes(Umatv(:,1:i),g(1:i,i)))*(yv-mtimes(Umatv(:,1:i),g(1:i,i))))/Nval;
        errorak(i) = error(i)*(1+2*i/Nval);
    end
        
end

function experiment2(a,b)
    Nest = 1000;
    Nval = 10000;
    u0 = randn(N,1);
    ny = normrnd(0,0.05,N,1);
    y0 = filter(b,a,u0);
    y = y0+ny;
end

function [Gzero] = G0(a,b,z)
    na = size(a,2);
    nb = size(b,2);
    suma=0.0;
    sumb=0.0;
    zk=1.0;
    for i=0:(na-1)
        suma= suma + a(i+1)*zk;
        zk = zk/z;
    end
    zk=1.0;
    for i=0:(nb-1)
        sumb = sumb + b(i+1)*zk;
        zk = zk/z;
    end
    Gzero = sumb / suma;
end

