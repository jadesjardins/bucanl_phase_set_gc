function [dN,R]=PSGCplugin_getPoints(Data,Ch1,M,varargin)

if nargin == 3
    Breaks=1;
elseif nargin == 4
    Breaks=varargin{1};
else
    error('error, wrong number of parameters')
end

% global variables
global w width c W RW deltaP T

% [b,a] are vectors which define a lowpass butter filter, 4th order
% The second parameter is a percentage of the nyqvist frequency it defines
% the cutoff point
% Filter k Hz, x*(T/2)=k, x=2*k/T
% Lowpass filter for complex demodulation algorithm
[b,a]=butter(4,2*width/T,'low');

dN=0;
R=0;
flag=0;

for B=1:length(Breaks)
    if B==length(Breaks)
        % X1 is channel of interest
        X1=double(Data(Ch1,Breaks(B):end));
    else
        X1=double(Data(Ch1, Breaks(B):(Breaks(B+1)-1)));
    end
    if length(X1) > (M*W+deltaP)
        %Partition the time series into K bins of length deltaP
        K=floor((length(X1)-M*W)/deltaP);
        interval=zeros(2,K);
        for i=1:K
            interval(1,i)= M*W+1+(i-1)*deltaP;
            interval(2,i)= M*W+1+i*deltaP;
        end
        % Initialize output variables
        R2=zeros(K,M);
        dN2=zeros(K,1);
        % index of sample
        t=1:length(X1);
        % windowing function
        window = hamming(length(X1));
        %Y - windowed sine (X) and cosine (Hilbert transform) components of X1, with frequency w shifted to 0
        Y1=window'.*X1.*sin(2*pi*w/T*t);
        Y2=window'.*X1.*cos(2*pi*w/T*t);
        
        % F is lowpass filter on Y to remove frequencies away from w
        F1=filtfilt(b,a,Y1);
        F2=filtfilt(b,a,Y2);
        % The power in the remaining signal is that from the band (w-k,
        % w+k) from the original signals
       

        % Calculate phase (P) of X1
        P=atan2(F2,F1);
        LP=length(P);
        P = unwrap(P);
        
        % Phase derivative
        dP=zeros(1,LP-2);
        for i=2:(LP-1)
            dP(i-1)=(P(i+1)-P(i-1))/2;
        end
        % Identify all the phase shifts
        % Count shift is number, N1 is the time index
        % If refract is greater than 0, the signal is in a refractory period and no
        % shifts are counted
        countShift=0;
        N1=0;
        refract=0;
        LdP = length(dP);
        for i=1:(LdP-1)
            if abs(dP(i)) < c && abs(dP(i+1)) > c && refract==0
                countShift=countShift+1;
                N1(countShift)=i;
                refract=RW;
            elseif refract>0
                refract=refract-1;
            end
        end
        
        % Assign shifts to dN
        Count1=zeros(1,K);
        for i=1:K
            Count1(i) = sum(N1 > M*W & N1 < interval(2,i));
        end
        dN2(1)=Count1(1);
        for i=2:K
            dN2(i)=Count1(i)-Count1(i-1);
        end
        % Find history R for this channel
        % Number of shifts in M previous windows of length W
        for i=1:K
            temp=zeros(2,M);
            for m=1:M
                temp(2,m)=interval(1,i)-(m-1)*W;
                temp(1,m)=interval(1,i)-m*W;        
            end
            for m=1:M
                R2(i,m) = sum(N1 > temp(1,m) & N1 < temp(2,m));
            end
        end
        if flag==0
            dN=dN2;
            R=R2;
            flag=1;
        else
            dN = [dN;dN2];
            R = [R;R2];
        end
    end
end
        
    
