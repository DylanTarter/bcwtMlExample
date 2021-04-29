function [header, stream] = BCWT(I,Qmin,Shift)
%% System Controls
wn = 'haar';
%% Pre-Processing
if(size(I,3) > 1) I = rgb2gray(I); end
[A,cH,cV,cD] = dwt2(I,wn);
C = round([zeros(size(A,1),size(A,2)),cH;cV,cD]);
C = (C<0).*(round(C-Shift)) +(C>0).*(round(C+Shift));
[A,cH,cV,cD] = dwt2(A,wn);
Cn = round([zeros(size(A,1),size(A,2)),cH;cV,cD]);
Cn = (Cn<0).*(round(Cn-Shift)) +(Cn>0).*(round(Cn+Shift));

Qb = floor(log2(abs(C)+(abs(C)<1))) + -1.*(abs(C)<1);
quad_max = @(block_struct) max(block_struct.data,[],'all');
MQD = blockproc(Qb, [2 2], quad_max);
Qb = Qb>=Qmin;

intvars = size(C,1)^2 + size(Cn,1)^2 + size(MQD,1)^2;
binvars = size(Qb,1)^2;
%disp(['Max Storage Cost: ', num2str(intvars),' Ints | ', num2str(binvars),' Bins = ', num2str(intvars+binvars)]);

%% Algorithm
%tic;
sqsize = size(MQD,2)/2;
stream = []; smax = 2;
levels = log2(size(I,2));
for l = 1:levels
    %disp(['Level: ',num2str(l)]);
    % Parse MQDn (not calculated)
    layers(l).MQD = MQD(1:smax*sqsize,1:smax*sqsize);
    layers(l).C = C(1:2*smax*sqsize,1:2*smax*sqsize);
    %layers(l).Qb = Qb(1:2*smax*sqsize,1:2*smax*sqsize);
    for x = 1:sqsize
        for y = 1:sqsize
            MQDt = [MQD(2*x-1,2*y-1),...
                    MQD(2*x-1,2*y);...
                    MQD(2*x,2*y-1),...
                    MQD(2*x,2*y)];
            Qg = max(max(MQDt));
            if(Qg>=Qmin)
            for i = 1:smax
                for j = 1:smax
                    if(MQDt(i,j) >= Qmin)
                        for n = 1:2
                            for m = 1:2
                                op = [(2*((2*x-1))-1+2*(i-1))+(n-1),...
                                    (2*((2*y-1))-1+2*(j-1))+(m-1)];
                                Ct = C(op(1),op(2));
                                if(Qb(op(1),op(2)))
                                    if(Ct >=0)
                                        stream = [1,stream];
                                    else
                                        stream = [0,stream];
                                    end
                                end
                                stream = [bn(B(abs(Ct)),MQDt(i,j),Qmin),stream];
                            end
                        end
                    end
                    stream = [bn(T(MQDt(i,j)),Qg,max(MQDt(i,j),Qmin)),stream];
                end
            end
            
            % Solving for MQDn
            Ct = [Cn(2*x-1,2*y-1),...
                  Cn(2*x-1,2*y);...
                  Cn(2*x,2*y-1),...
                  Cn(2*x,2*y)];
            Qt = floor(log2(abs(Ct)+(abs(Ct)<1))) + -1.*(abs(Ct)<1);
            MQD(x,y) = max(max(max(Qt)),Qg);
            stream = [bn(T(Qg),MQD(x,y),max(Qg,Qmin)),stream];
            Qb(2*x-1:2*x,2*y-1:2*y) = (Qt>=Qmin);
            C(2*x-1:2*x,2*y-1:2*y) = Ct;
            end
        end
    end
    sqsize = sqsize / 2;
    if(sqsize < 1)
        sqsize = 1;
        smax = 1;
    else
        [A,cH,cV,cD] = dwt2(A,wn);
    end
    Cn = round([zeros(size(A,1),size(A,2)),cH;cV,cD]);
    Cn = (Cn<0).*(round(Cn-Shift)) +(Cn>0).*(round(Cn+Shift));
end

% Final Level Handler
%toc;
%% Create the Header
qmax = Qg;
header = [];
header = [header,dec2bin(Qmin,4) - '0'];
header = [header,dec2bin(qmax,4) - '0'];
header = [header,dec2bin(levels,4) - '0'];
% Final Level (For single point Compression)
header = [header,dec2bin(A,16) - '0'];

end
%%
function bin = B(dec) 
    bin = de2bi(dec,16);
    if(length(bin) < 16) 
        bin = [bin,zeros(1,16-length(bin))];
    end
    bin = fliplr(bin);
end
function out = bn(bin,m,n)
    out = bin(16-m:16-n);
end
function bin = T(pos) 
    bin = zeros(1,16);
    bin(16-pos) = 1;
end