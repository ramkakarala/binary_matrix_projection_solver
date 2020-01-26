

% create 3D matrix
N = 4;
M = zeros(N,N,N);

Z = [ 1,1,1,1 ; 0,0,1,0 ; 0,1,0,0 ; 1,1,1,1 ];
R = [ 1,1,1,1 ; 1,0,0,1 ; 1,1,1,0 ; 1,0,0,1 ];
FIT = [ 1,0,0,1 ; 1,0,0,0 ; 1,0,0,0 ; 1,0,0,0 ];
    
M(1,:,:) = Z; 
M(:,1,:) = R;
M(:,:,1) = FIT;



figure;
subplot(2,2,1)
imagesc(squeeze(sum(M,1)>0)) ;
colormap(gray)
title('Sum over 1st');

subplot(2,2,2)
imagesc(squeeze(sum(M,2)>0)) ;
colormap(gray)
title('Sum over 2nd');

subplot(2,2,3)
imagesc(squeeze(sum(M,3)>0)) ;
colormap(gray)
title('Sum over 3rd');


% https://www.mathworks.com/help/matlab/ref/sphere.html


[x,y,z] = sphere ;
x = 0.5*x;
y = 0.5*y;
z = 0.5*z;
figure;
for ii=1:N,
    for jj=1:N,
        for kk=1:N,
            if (M(ii,jj,kk)==1),
                surf(x+ii,y+jj,z+kk) % centered at (3,-2,0) 
                hold on
            end
        end
    end
end

% https://www.mathworks.com/help/matlab/ref/rotate3d.html
% MANY exampels there of rotating only parts of a plot, or using buttons
h = rotate3d;
%h.RotateStyle = 'box';
h.Enable = 'on';



%