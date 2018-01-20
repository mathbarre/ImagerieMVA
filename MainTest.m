%im = (imread('TP/im/scotland_house.png')); 
%pathImA = 'TP/im/flower3.jpeg';
%pathImB = 'TP/im/flower2.bmp';
%pathImA = 'TP/im/forest_summer.jpg';
%pathImB = 'TP/im/forest_autumn.jpg';
%pathImB = 'TP/im/dicaprio.jpg';
%pathImA = 'TP/im/avatar.jpg';
%pathImA = 'TP/im/ice_dark.jpg';
%pathImB = 'TP/im/ice_clear.jpg';
%pathImB = 'TP/im/sun_ice.jpg';
%pathImA = 'TP/im/sunset_ice.jpg';
%pathImB = 'TP/im/flower1.bmp';
%pathImA = 'TP/im/flower2.bmp';
%pathImA = 'TP/im/ice_dark.jpg';
%pathImB = 'TP/im/sunset2.jpg';
%pathImA = 'TP/im/avatar.jpg';
%pathImB = 'TP/im/femme.jpg';
%pathImA = 'TP/im/adachi_moins_vert.jpeg';
%pathImB = 'TP/im/adachi_vert.jpg';
pathImA = 'TP/im/scotland_house.png';
pathImB = 'TP/im/scotland_plain.png';
%pathImB ='TP/im/scotland_house.png';
%pathImA = 'TP/im/scotland_house.png';
%pathImB = 'TP/im/tropique.jpg'

%pathImA = 'TP/im/fuji_automne.jpg';
%pathImB = 'TP/im/fuji_hiver.jpg';

%pathImA = 'Results/scotland_plain_on_scotland_house_eps_3_lambda_0.001_R_60_regrain.png';
nameA = strsplit(pathImA,{'/','.'});
nameA = nameA{3};
nameB = strsplit(pathImB,{'/','.'});
nameB = nameB{3};
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

