function [R2M,L2M,T2M,Hs,RM]=croscor_20180904_modv2(A2,dpha,dt,dc)

% Copyright 2017 Rafael Almar (IRD, France)- rafael.almar@ird.fr
% A2: input matrix (spatio-temporal matrix) [ntim ncol]
% dpha: (For method 1) fixed time shift
% dc: (For method 1) maximum spatial shift
% dt: time step
% res: spatial resolution of the calculation
% nmax (For methods 2 3 4) maximum time shift
% methodcor :(For method 1) choice of peak (max=1, median=2)
% method: 1 = calculate using time, 2 = calculate time shift (fast method using fft)
%          3 = calculate time shift (slower method using time)
%          4 = calculate between two vectors only, calculate time shift
% ecart: pixel distance that characterizes the spread of the correlation peak
%        (associated with a value (ex: val = 2/3 of the peak value)
% pval: the value of the p-value. Each p-value is the probability of getting a
%       correlation as large as the observed value by random chance, when the
%       true correlation is zero. If P(i,j) is small, say less than 0.05,
%       then the correlation R(i,j) is significant. (associated with a
%       significance level of 0.05 or 0.01)
% SigExpl: proportion of the total signal explained by the linear regression model
A2 = double(A2(:,:,1));
% A2 = detrend(double(abs(diff(A2))));
A2 = detrend((A2));
dc = round(dc/2).*2;
res = 1; %Resolution for cross-correlation

[nt,nc,cc]=size(A2);

n = floor(dpha/dt);
R2(1:dc-1) = 0;

for ic=1+dc/2:res:nc-dc/2
    for lc=1:dc-1
        % variable in time
        R7 = corrcoef(A2(n+1:nt,ic-dc/2),A2(1:(nt-(n)),ic-dc/2+lc-1));
        R2(lc)=R7(2,1);
    end 
    R2M(ic-dc/2,1:dc-1) = R2;
end % ic

%Wavelength
for i=1:size(R2M,1)
    [hs,htiers,hrms,trms,Hmax,h,Tp,sl,tcrete,t]=Wave_Char(R2M(i,:),1,0,2);
    L2M(i) = trms ;
end

%Wave period
for  i=1+dc/2:nc-dc/2
    try
        [hs,htiers,hrms,trms,Hmax,h,Tp,sl,tcrete,t]=Wave_Char(detrend(A2(1:min([500 round(length(A2(:,1)))]),i)),1,0,2);
        T2M(i) = trms.*dt;
        Hs(i) = hs;
    catch
        T2M(i)= NaN;
        Hs(i) = NaN;
    end
end
RM = R2M;
[lm R2M]=max(R2M(:,:)');R2M = R2M./dpha;
end
