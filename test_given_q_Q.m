function [  ] = given_q_Q_test_all(  )
%  given image compression parameters q and Q, test compression efficiency
%   
addpath(genpath('hdrvdp-2.2.1'));
% hdr_folder = 'D:\\work\\data\\subset\\';
hdr_folder = 'codec_result\\JPEG-XT-HDR20\\';
hdr_relative_folder ='JPEG-XT-HDR20\\';
hdrFileList = dir([hdr_folder '*.hdr']);
%% print all files
fnameList = cell(length(hdrFileList), 1);
for fno=1:1:length(hdrFileList)
    [~,fname,~] = fileparts(hdrFileList(fno).name);
    fnameList{fno} = fname;
end

epsilon = 1000;
% q = 90; Q = 90;

codec = cell(3,1);
codec{1} = 'jpghdr_raw.exe';
codec{2} = 'jpghdr_no_decomposition.exe';
codec{3} = 'jpghdr_with_decomposition.exe';

% tmos = { 'Reinhard2002'};
% tmos = {'Reinhard2002', 'FWLS', 'GIF'};

tmos = {'mantuik06'};

q = 65;
paras_pair = [q,  10;...
              q,  35; ...
              q,  65; ...
              q,  90];
% paras_pair = [q, 90];
qQres = zeros(length(hdrFileList),12);
qQfile = zeros(length(hdrFileList),12);
for k=1:1:size(tmos,2)
    result = zeros(length(hdrFileList), 8);
    fprintf('TMO: %s\n', tmos{k});
    for fno = 1:1:length(hdrFileList)  
        [~,fname,~] = fileparts(hdrFileList(fno).name);
        hdr_fullpath = sprintf('%s%s.hdr', hdr_folder, fname);
        hdr_relativepath = sprintf('%s%s.hdr', hdr_relative_folder, fname);
        hdr = hdrread(hdr_fullpath);
        %% compress ldr with the given q 
        cd('.\\codec_result');
        jpg_comp_file = sprintf('%s_q%d.jpg', tmos{k}, q);
        comp_cmd = sprintf('standard_jpeg.exe -q %d .\\%s\\%s.pbm %s', ...
                            q, tmos{k}, fname, jpg_comp_file);
        system(comp_cmd);
        ldr = imread(jpg_comp_file);
        cd('..');

        ldr = double(ldr)/255;
        ldr = ldr .* ldr; % inverse gamma correction

        hy = rgb2Y(hdr); ly = rgb2Y(ldr);
        pos_hy = remove_zero_luminance(hy);
        pos_ly = remove_zero_luminance(ly);

        loghy = log2(pos_hy); logly = log2(pos_ly); % use log2, not log10
        
        for qidx=1:1:size(paras_pair,1)
             Q = paras_pair(qidx, 2);
             
           %% use default JPEG_HDR
            fprintf('###### jpeg_raw   ######\n');
            cd('.\\codec_result');
            ldr_pbm = sprintf('.\\%s\\%s.pbm', tmos{k}, fname);
            jpg_out_file = 'output.jpg';
            hdr_out_file = 'hdr_recovery.hdr';
            % jpghdr_raw.exe
            enc_cmd = sprintf('%s -t %s -q %d -Q %d %s %s',codec{1}, ldr_pbm, ...
                            q, Q, hdr_relativepath, jpg_out_file);
            system(enc_cmd);
            tmp = dir(jpg_out_file);
            result(fno, 1) = tmp.bytes/1024;
            dec_cmd = sprintf('%s %s %s', codec{1}, jpg_out_file, hdr_out_file);
            system(dec_cmd);
            hdr_recovery = hdrread(hdr_out_file);
            cd('..');
            res = hdrvdp( hdr, hdr_recovery, 'rgb-bt.709', 30, { 'surround_l', 10^-5 } );
            result(fno, 2) = res.Q;
            
            fprintf('#########################\n');
            qQfile(fno, 3*qidx - 2) = tmp.bytes/1024;
            qQres(fno, 3*qidx - 2) = res.Q;
            
            % without decomposition
            [loghb, loglb] = generate_base_layers(loghy, logly, 0);
            [itmf, min_loglb, max_loglb] = estimate_itmf(loghb, loglb, epsilon);
            fprintf('###### jpeg_no_decomposition   ######\n');
            cd('.\\codec_result');
            ldr_pbm = sprintf('.\\%s\\%s.pbm', tmos{k}, fname);
            jpg_out_file = 'output.jpg';
            hdr_out_file = 'hdr_recovery.hdr';
%             %jpghdr_ours.exe
            enc_cmd = sprintf('%s -t %s -q %d -Q %d %s %s',codec{2}, ldr_pbm, ...
                            q, Q, hdr_relativepath, jpg_out_file);
            system(enc_cmd);
            tmp = dir(jpg_out_file);
            result(fno, 3) = tmp.bytes/1024;
            dec_cmd = sprintf('%s %s %s', codec{2}, jpg_out_file, hdr_out_file);
            system(dec_cmd);
            hdr_recovery = hdrread(hdr_out_file);
            cd('..');
            res = hdrvdp( hdr, hdr_recovery, 'rgb-bt.709', 30, { 'surround_l', 10^-5 } );
            result(fno, 4) = res.Q;
            
            fprintf('#########################\n');
            qQfile(fno, 3*qidx -1) = tmp.bytes/1024;
            qQres(fno, 3*qidx - 1) = res.Q;
            
            if strcmp(tmos{k}, 'Reinhard2002')
            else
                %% with decomposition
                [loghb, loglb] = generate_base_layers(loghy, logly, 1);
                [itmf, min_loglb, max_loglb] = estimate_itmf(loghb, loglb, epsilon);                
                
                fprintf('###### jpeg_with_decomposition   ######\n');
                cd('.\\codec_result');
                ldr_pbm = sprintf('.\\%s\\%s.pbm', tmos{k}, fname);
                jpg_out_file = 'output.jpg';
                hdr_out_file = 'hdr_recovery.hdr';
                % jpghdr_ours.exe
                enc_cmd = sprintf('%s -t %s -q %d -Q %d %s %s',codec{3}, ldr_pbm, ...
                                q, Q, hdr_relativepath, jpg_out_file);
                system(enc_cmd);
                tmp = dir(jpg_out_file);
                result(fno, 5) = tmp.bytes/1024;
                dec_cmd = sprintf('%s %s %s', codec{3}, jpg_out_file, hdr_out_file);
                system(dec_cmd);
                hdr_recovery = hdrread(hdr_out_file);
                cd('..');
                res = hdrvdp( hdr, hdr_recovery, 'rgb-bt.709', 30, { 'surround_l', 10^-5 } );
                result(fno, 6) = res.Q;
                
                
                qQres(fno, 3*qidx) = res.Q;
                qQfile(fno, 3*qidx) = tmp.bytes/1024;
            end
        end
%     save(sprintf('result_tmo_%d.mat', k), 'result');
    end
end
disp('Complete.\n');
end





function [loghb, loglb] = generate_base_layers(loghy, logly, with_decomposition)
    if with_decomposition
        loglb = fast_2D_smoother(logly, logly, 1, 3);
%         loghb = fast_2D_smoother(loghy, loghy, 1, 3);
        loghb = loghy - (logly - loglb);
    else
        loglb = logly;
        loghb = loghy;
    end
end

% estimate inverse tone mapping curve using the logarithm of base layers
function [itmf, min_loglb, max_loglb] = estimate_itmf(loghb, loglb, epsilon, min_loglb, max_loglb)
  n = size(loghb,1)*size(loghb, 2); % number of pixels

% split interval [quantitized_loglb_min quantitized_loglb_max] into k intervals
    k = 256; % k = 256, 512, 1024
    if nargin == 3
        % normalize and quantitize
        min_loglb = min(loglb(:));
        max_loglb = max(loglb(:));

%         min_loglb = 1.025*min_loglb;
%         max_loglb = 0.975*max_loglb;
    elseif nargin == 5
        
    end
    int_range = (max_loglb - min_loglb) / (k-1);
    loglb(loglb < min_loglb) = min_loglb;
    loglb(loglb > max_loglb) = max_loglb;
    quantitized_loglb = floor((loglb - min_loglb)/int_range + 0.5);

    quantitized_loglb = uint16(quantitized_loglb); % if k >256, it should be uint16

    quantitized_loglb = quantitized_loglb(:);
    k = uint16(k);
% compute weight LUT
    w = zeros(k, 1);
    for i=1:1:(k/2) % 128
        w(i) = (i-1);
%         w(i) = 1;
    end
    for i=(k/2+1):1:k
        w(i) = (k - i);
%         w(i) = 1;
    end

    Y = loghb(:);

    for i = 1:1:n
        Y(i) = Y(i)*w(uint16(quantitized_loglb(i))+1);
    end

    M = zeros(k-2, k);
    for i=1:1:(k-2)
        M(i, i) = 1*sqrt(w(i+1)); 
        M(i, i+1) = -2*sqrt(w(i+1));
        M(i, i+2) = 1*sqrt(w(i+1));
    end

% H'*H is a diagonal matrix
% H'*W*H = diag(w0*n0, w1*n1, ..., w255*n255)
    D = zeros(k,k);
    HY = zeros(k,1);

    for i=1:1:k
        temp = find(quantitized_loglb == (i-1));
        HY(i) = sum(Y(temp));
        n_i = length(temp);
%         fprintf('%d\n', n_i);
        D(i,i) = w(i)*n_i;
    end

    A = D + epsilon* (M')*M;

    itmf = A \ (HY);
    
    %% write itmf to file
    fid = fopen('.\\codec_result\\itmf.txt', 'w');
    fprintf(fid, '%.8f\n%.8f\n', min_loglb, max_loglb);
    for k=1:1:256 % 256 bins
        fprintf(fid, '%.8f\n', itmf(k));
    end
    fclose(fid);
end

% the luminance of HDR image is always larger than 0, force 0 luminance to
% the smallest positive luminance to reduce estimation error
function [positive_lum] = remove_zero_luminance(linear_lum)
    smallest_positive_lum = 65536;
    for i=1:1:size(linear_lum,1)
        for j=1:1:size(linear_lum,2)
            if (linear_lum(i,j)<smallest_positive_lum) && (linear_lum(i,j)>0)
                smallest_positive_lum = linear_lum(i,j);
            end
        end
    end
    positive_lum = linear_lum;
    for i=1:1:size(linear_lum,1)
        for j=1:1:size(linear_lum,2)
            if (positive_lum(i,j)<=0)
                positive_lum(i,j) = smallest_positive_lum;
            end
        end
    end
end

function [Y] = rgb2Y(I)
    Y = 0.299*I(:,:,1) + 0.587*I(:,:,2) + 0.114*I(:,:,3);
end