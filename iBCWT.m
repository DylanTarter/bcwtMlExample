%% unpacking the header1
function I2 = iBCWT(header,stream)
header1 = header;
qmin = bin2dec(num2str(header1(1:4)));header1(1:4) = [];
qmax = bin2dec(num2str(header1(1:4)));header1(1:4) = [];
mlevels = bin2dec(num2str(header1(1:4)));header1(1:4) = [];
A = bin2dec(num2str(header1(1:16)));header1(1:16) = [];

%disp(['qmin = ',num2str(qmin),' | qmax = ',num2str(qmax)]);
%disp(['MQD Levels = ',num2str(mlevels)]);

% MQD - QL = bit stuff -> QL = MQD - bit stuff
clc;
% added a 0 at end so the loops dont break on the case of the last thing is
% nothing
stream1 = [stream,0];
MQD = qmax;
level = 1;
layer(1).MQD = MQD;

% Pull off the first guys
MQD_n = zeros(1,1);
C = zeros(2,2);
QL = qmax;
QL0 = QL;
piece = [];
while(stream1(1) ~= 1 && QL >= qmin)
    piece = [piece, stream1(1)]; stream1(1) = [];
    QL = QL - 1;
end
if(isequal(piece,zeros(1,QL0-qmin+1))) QL = 0; else stream1(1) = [];end
QL_t = QL;
while(stream1(1) ~= 1 && QL_t >= qmin)
    piece = [piece, stream1(1)]; stream1(1) = [];
    QL_t = QL_t - 1;
end
if(isequal(piece,zeros(1,QL-qmin+1))) QL_t = 0; else stream1(1) = [];end
BPS = (QL_t-qmin)+1;
for c = 1:2
    for v = 1:2
        piece = [];
        for k = 1:(QL_t-qmin)+1
            piece = [piece, stream1(1)];
            stream1(1) = [];
        end
        str = num2str([piece,zeros(1,qmin)]);
        str(isspace(str)) = '';
        num = bin2dec(str);
        if(num~=0)
            my_sign = stream1(1);
            stream1(1) = [];
        if(my_sign==0) num = num * -1; end
        end
        C(3-c,3-v) = num;
    end
end
MQD = QL;
layer(1).MQD = MQD;
layer(1).C = C;
% pop off layer 1
cV = C(floor(end/2)+1:end,1:floor(end/2));
cH = C(1:floor(end/2),floor(end/2)+1:end);
cD = C(floor(end/2)+1:end,floor(end/2)+1:end);
A = idwt2(A,cH,cV,cD,'haar');

%%
sq_size_p = 2^(level-1);
sq_size = 2^(level);
%tic;
% x , y = analyze QL's
% i , y = create MQD next
for l = 1:mlevels-1
    %disp(['Level: ',num2str(l)]);
    MQD_n = zeros(sq_size,sq_size);
    C = zeros(sq_size*2,sq_size*2);
    xcoord = 1:2:sq_size;
    ycoord = 1:2:sq_size;
    cxcoord = 1:2:sq_size*2;
    cycoord = 1:2:sq_size*2;
    for x = 1:sq_size_p
        for y = 1:sq_size_p
            QL = MQD(end-x+1,end-y+1);
            QL0 = QL;
            MQD_n1 = zeros(2,2);
            if(QL ~= 0)
            piece = [];
            while(stream1(1) ~= 1 && QL >= qmin)
                piece = [piece, stream1(1)]; stream1(1) = [];
                QL = QL - 1;
            end
            if(isequal(piece,zeros(1,QL0-qmin+1))) QL = 0; else stream1(1) = []; end
            C_unit = zeros(4,4);
            cux = 1:2:4;
            cuy= 1:2:4;
            if(QL~=0)
            for i = 1:2
                for j = 1:2
                    QL_t = QL;
                    piece = [];
                    while(stream1(1) ~= 1 && QL_t >= qmin)
                        piece = [piece, stream1(1)]; stream1(1) = [];
                        QL_t = QL_t - 1;
                    end
                    if(isequal(piece,zeros(1,QL-qmin+1))) QL_t = 0; else stream1(1) = [];end
                    C_block = zeros(2,2);
                    if(QL_t >= qmin)
                        BPS = (QL_t-qmin)+1;
                        for c = 1:2
                            for v = 1:2
                                piece = [];
                                for k = 1:BPS
                                    piece = [piece, stream1(1)];
                                    stream1(1) = [];
                                end
                                str = num2str([piece,zeros(1,qmin)]);
                                str(isspace(str)) = '';
                                num = bin2dec(str);
                                if(num~=0)
                                    my_sign = stream1(1);
                                    stream1(1) = [];
                                    if(my_sign==0) num = num * -1; end
                                end
                                C_block(3-c,3-v) = num;
                            end
                        end
                        MQD_n1(end-i+1,end-j+1) = QL_t;
                    else
                        MQD_n1(end-i+1,end-j+1) = 0;
                    end
                    C_unit(cux(3-i):cux(3-i)+1,cuy(3-j):cuy(3-j)+1) = C_block;
                end
            end
            end
                C(cxcoord(end-x*2+1):cxcoord(end-x*2+1)+3,cycoord(end-y*2+1):cycoord(end-y*2+1)+3) = C_unit;
            else
                C(cxcoord(end-x*2+1):cxcoord(end-x*2+1)+3,cycoord(end-y*2+1):cycoord(end-y*2+1)+3) = zeros(4,4);
            end
            MQD_n(xcoord(end-x+1):xcoord(end-x+1)+1,ycoord(end-y+1):ycoord(end-y+1)+1) = MQD_n1;
        end
    end
    layer(l+1).C = C;
    layer(l+1).MQD = MQD_n;
    sq_size_p = sq_size;
    sq_size = 2^(l+1);
    MQD = MQD_n;
    
    % Inverse Transform
    cV = C(floor(end/2)+1:end,1:floor(end/2));
    cH = C(1:floor(end/2),floor(end/2)+1:end);
    cD = C(floor(end/2)+1:end,floor(end/2)+1:end);
    A = idwt2(A,cH,cV,cD,'haar');
end
%toc;

%% Display
% subplot(1,2,1);
% I1 = I;
% imagesc(I1);
% kB1 = imfinfo(filename).FileSize/1000;
% xlabel(['Bit Count: ',num2str(kB1),'kB']);
% subplot(1,2,2);
 I2 = uint8(A);
% imagesc(I2);
% kB2 = (size(stream,2)+size(header,2))/8000;
% xlabel(['Bit Count: ',num2str(kB2),'kB']);
% colormap('gray');
% sgtitle(['BCWT Compression for Qmin=',num2str(qmin),' | Compression: '...
%     ,num2str(size(stream,2)/(size(I2,1)*size(I2,1)*8)),'bpp | Max Diff: ',num2str(max(max(abs(I1-I2))))]);
% exp(end+1) = psnr(I2,I1);
% exp1(end+1) = kB1/kB2;