%Experiment1
figure('Name','Experiment 1')
experiment1(1000,5000,1000,0.1,[0.999 0.95 0.6],0.1,0.1,1.0,1);

%Experiment2
figure('Name','Experiment 2')
experiment2(1000,5000,1000,0.1,0.6,0.1,0.1,1.0,[1 2 5]);

%Experiment3
figure('Name','Experiment 3')
experiment3(1000,5000,1000,0.01,0.001,1.0);
hi0 =0;
hni=0;

function experiment1(tests,N,R0,fgen,fnoise,sigmai0,sigmani,sigmanu,IVstride)
    LS = zeros(tests,3);
    IV = zeros(tests,3);
   
    for k=1:3
        for  j=1:tests
            [i,u,i0,ni,hi0,hni] = generatefilter(R0,fgen, fnoise(k), sigmai0, sigmani,sigmanu,N);%TODO fix outputs
            LS(j,k)=LSstimator(i,u);
            IV(j,k)=IVstimator(i,u,IVstride);
           
        end
        subplot(3,3,3+k);
        [r,lags] = xcorr(i0,'unbiased');
        plot(lags,r,'-o');
        hold on
        xlim([-0.1 10])
        [r,lags] = xcorr(ni,'unbiased');
        plot(lags,r,'-o');
        xlabel('Lag');
        ylabel('Auto-correlation');
        legend('i0','ni');
        hold off
        subplot(3,3,6+k);
        plot(hi0)
        hold on
        plot(hni)
        xlabel('f/fs');
        ylabel('Filter[dB]');
        xlim([0 0.5])
        legend('i0','ni');
        hold off
    end
    subplot(3,3,[1 2 3])
    hold on
    histogram(LS(:,1),'FaceColor',[0,0,0],'normalization','pdf')%BLACK
    for k=1:3
        histogram(IV(:,k),'FaceColor',[k==1,k==2,k==3],'normalization','pdf')
    end
    ylim([0 0.05]);
    ylabel('pdf(R)')
    xlim([480 1500]);
    legend('LS','IV(.999)','IV(.95)','IV(.6)');
    hold off
end

function [IV]=experiment2(tests,N,R0,fgen,fnoise,sigmai0,sigmani,sigmanu,IVstride)
    LS = zeros(tests,1);
    IV = zeros(tests,3);
    for  j=1:tests
        [i,u,i0,ni,hi0,hni] = generatefilter(R0,fgen, fnoise, sigmai0, sigmani,sigmanu,N); %TODO fix outputs
        LS(j)=LSstimator(i,u);
        IV(j,1)=IVstimator(i,u,IVstride(1));
        IV(j,2)=IVstimator(i,u,IVstride(2));
        IV(j,3)=IVstimator(i,u,IVstride(3));
    end
    for k=1:3
        subplot(3,3,3+k);
        [r,lags] = xcorr(i0,'unbiased');
        plot(lags,r,'-o');
        hold on
        xlim([-0.1 10])
        [r,lags] = xcorr(ni,'unbiased');
        plot(lags,r,'-o');
        xlabel('Lag');
        ylabel('Auto-correlation');
        legend('i0','ni');
        hold off
        subplot(3,3,6+k);
        plot(hi0)
        hold on
        plot(hni)
        xlabel('f/fs');
        ylabel('Filter[dB]');
        xlim([0 0.5]);
        legend('i0','ni');
        hold off
    end
    subplot(3,3,[1 2 3])
    histogram(LS,'FaceColor',[0,0,0],'normalization','pdf')%BLACK
    hold on
    histogram(IV(:,1),'FaceColor',[1,0,0],'normalization','pdf')%RED
    histogram(IV(:,2),'FaceColor',[0,1,0],'normalization','pdf')%GREEN
    histogram(IV(:,3),'FaceColor',[0,0,1],'normalization','pdf')%BLUE
    ylim([0 0.05]);
    ylabel('pdf(R)')
    xlim([480 1500]);
    legend('LS','IV(1)','IV(2)','IV(5)');
    hold off
end

function experiment3(tests,N,R0,sigmai0,sigmani,sigmanu)
    LS = zeros(tests,1);
    EIV = zeros(tests,1);
    for  j=1:tests
        [i,u] = generate(R0,sigmai0, sigmani,sigmanu,N);%TODO fix outputs
        LS(j) = LSstimator(i,u);
        EIV(j) = EIVstimator(i,u,sigmani,sigmanu);
    end
    histogram(LS,'FaceColor',[0,1,0],'Normalization','pdf')%GREEN
    hold on
    histogram(EIV,'FaceColor',[1,0,0],'Normalization','pdf')%RED
    legend('LS','EIV');
    hold off
end

function [R] = LSstimator(i,u)
    R =(transpose(i)*u)/(transpose(i)*i);
end

function [R] = IVstimator(i,u,s)
    i2 = transpose(circshift(i,s));
    R =(i2*u)/(i2*i);
end

function [R] = EIVstimator(i,u,sigmani,sigmanu)
    su = sigmanu^2;
    si = sigmani^2;
    i2 = ((transpose(u)*u)/su)-((transpose(i)*i)/si);
    iu = transpose(i)*u;
    R =(i2+sqrt((i2^2)+(4*(iu^2)/(si*su))))/(2*iu/su);
end

%TODO fix outputs
function [i,u,i0,ni,hi0,hni] = generatefilter(R0,fgen, fnoise, sigmai0, sigmani,sigmanu,N)
    e1=randn(N,1);%current source
    [bgen,agen] = butter(1,fgen); %generates butterworth first order filter
    hi0 = freqz(bgen,agen);
    i0 = filter(bgen,agen,e1); %applies the filter to the current
    i0=i0/(std(i0)/sigmai0);
    %generating current noise
    e2=randn(N,1);
    [bnoise,anoise] = butter(2,fnoise);%second order butterwoth filter
    hni = freqz(bnoise,anoise);
    ni = filter(bnoise,anoise,e2);
    nu = normrnd(0,sigmanu,N,1);
    ni=ni/(std(ni)/sigmani);
    nu=nu/(std(nu)/sigmanu);
    u0 = i0 * R0;
    i=i0+ni;
    u=u0+nu;
end

%TODO fix outputs
function [i,u] = generate(R0, sigmai0, sigmani,sigmanu,N)
    i0 = normrnd(0,sigmai0,N,1);
    %generating current noise
    ni = normrnd(0,sigmani,N,1);
    nu = normrnd(0,sigmanu,N,1);
    u0 = i0 * R0;
    i=i0+ni;
    u=u0+nu;
end
