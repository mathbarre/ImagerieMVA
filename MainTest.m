%im = (imread('im/scotland_house.png')); 
%pathImA = 'im/flower3.jpeg';
%pathImB = 'im/flower2.bmp';
%pathImA = 'im/forest_summer.jpg';
%pathImB = 'im/forest_autumn.jpg';
%pathImB = 'im/dicaprio.jpg';
%pathImA = 'im/avatar.jpg';
%pathImA = 'im/ice_dark.jpg';
%pathImB = 'im/ice_clear.jpg';
%pathImB = 'im/sun_ice.jpg';
%pathImA = 'im/sunset_ice.jpg';
%pathImB = 'im/flower1.bmp';
%pathImA = 'im/flower2.bmp';
%pathImA = 'im/ice_dark.jpg';
%pathImB = 'im/sunset2.jpg';
%pathImA = 'im/avatar.jpg';
%pathImB = 'im/femme.jpg';
%pathImA = 'im/adachi_moins_vert.jpeg';
%pathImB = 'im/adachi_vert.jpg';
pathImA = 'im/scotland_house.png';
pathImB = 'im/scotland_plain.png';
%pathImB ='im/scotland_house.png';
%pathImA = 'im/scotland_house.png';
%pathImB = 'im/tropique.jpg'

%pathImA = 'im/fuji_automne.jpg';
%pathImB = 'im/fuji_hiver.jpg';

%pathImA = 'Results/scotland_plain_on_scotland_house_eps_3_lambda_0.001_R_60_regrain.png';
nameA = strsplit(pathImA,{'/','.'});
nameA = nameA{2};
nameB = strsplit(pathImB,{'/','.'});
nameB = nameB{2};
im = (imread(pathImA));
imB = (imread(pathImB)); 

[numRows,numCols,~] = size(im);


nbSuperPixelsWantedA = round(numRows*numCols/500);
%nbSuperPixelsWantedB = round(numRowsB*numColsB/500)+10;
nbSuperPixelsWantedB = nbSuperPixelsWantedA+20;

lambda = 0.001;
epsilon = 3;
R = 50;

% possible colorspace arg : 'rgb','opp','lms','lab','hsv'
resSCT = SCT(im,imB,epsilon,lambda,R,0.7,15,nbSuperPixelsWantedA,nbSuperPixelsWantedB,10,0.1,'rgb');

figure;imshow(resSCT,[]);
Final = regrain(double(im)/255,resSCT);
figure;imshow(Final,[]);

%imwrite(TransformedImage/255, strcat('Results/',nameB,'_on_',nameA,'_eps_',num2str(epsilon),'_lambda_',num2str(lambda),'_R_',num2str(R),'.png'));
%imwrite(Final, strcat('Results/',nameB,'_on_',nameA,'_eps_',num2str(epsilon),'_lambda_',num2str(lambda),'_R_',num2str(R),'_regrain.png'));

