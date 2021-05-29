function result = colourMatrix(filename)                         
                                % Image Transformation %
                                
baseline=imread('org_1.png'); % Reads the image, Here the org_1 image was chosen as our baseline image.
%figure,imshow(baseline),title('Baseline')
baselinebw  =~ im2bw(baseline,0.5); % Converts the image to binary (black and white).
erobaselinebw = imerode(baselinebw,ones(7)); % Erosion performed for denoise purpose.
dilbaselinebw = imdilate(erobaselinebw,ones(7)); % Dilation was used to sharp the eroded image.
%figure,imshow(dilbaselinebw),title('Binary Baseline with filter applied').

                              % Find Circle Regions %
           
CC = bwconncomp(dilbaselinebw); % Finds the connected components of the image.
z = regionprops(dilbaselinebw,'Area','Centroid'); % Calculates stats and find useful insights (that will be used as thresholds to help the circle recognition).
Areas = [z.Area]; % Making an array saving the area mass of each object.
%Areas
                                % Isolate Circles %
                    
isolatedCircles=zeros(CC.ImageSize); % Making a zero array, with equal size as the image.
for p=1:CC.NumObjects  %Loop through each object, identified from bwconncomp.
    if z(p).Area < max(Areas) % Condition statement that returns only objects that have smaller area mass than the maximum area.
        isolatedCircles(CC.PixelIdxList{p}) = 1; %Set the image.
    end

end
%figure,imshow(isolatedCircles,[0 1]) %Display with fixed intestiy range.
%title("Baseline's Isolated Circles")
labeledCircles = bwlabel(isolatedCircles, 8); %Labels each blob in the image, in our case, circles.
CircleMeasurements = regionprops(labeledCircles, 'all'); %Calculates stats of the circles.
allCircleCentroids = [CircleMeasurements.Centroid]; %Saving the four centroids of the four circles, in an array.
FixedPoints=reshape(allCircleCentroids,2,[]); %Sets the centroids as moving points.
tFixed = transpose(FixedPoints); %Transposes the matrix so it will have the form [xi yi; xi+1 yi+1; and so on].

                        % Load in the input image to be transformed %
       
inputImg=imread(filename) ;  
%figure,imshow(inputImg),title('Input Image')
inputBW =~ im2bw(inputImg,0.5); 
inputBWfill=imfill(inputBW,'holes'); % Denoises the image, filling noise holes.
inputMedfilt = medfilt2(inputBWfill); % Removing salt and pepper noise using the non-linear median filter.
%figure,imshow(inputMedfilt), title('Binary Input Image with filter applied')
                               % Isolate Circles %
                                 
CC2 = bwconncomp(inputMedfilt);
s = regionprops(inputMedfilt,'Area','Centroid');
Areas2 = [s.Area];
isolatedCircles2 = zeros(CC2.ImageSize);
for p = 1:CC2.NumObjects 
    if s(p).Area < max(Areas2)
        isolatedCircles2(CC2.PixelIdxList{p}) = 1; 
    end
end
%figure,imshow(isolatedCircles2,[0,1]), title("Input's Image Isolated Circles")

labeledCircles2 = bwlabel(isolatedCircles2, 8);
CircleMeasurements2 = regionprops(labeledCircles2,'Area','Centroid');
allCircleCentroids2 = [CircleMeasurements2.Centroid];
MovingPoints=reshape(allCircleCentroids2,2,[]); %Set Centroids as moving points.
tMoving = transpose(MovingPoints);

                         % Applying the geometric transformation. %

tform=fitgeotrans(tMoving,tFixed, 'Projective'); 

                                 % Size Correction %

R=imref2d(size(baseline)); %Resizes the transformed image, with reference the baseline.
transformedImg = imwarp(inputImg,tform,'OutputView',R);
%figure,imshow(transformedImg),title('Transformed Image')

                            % Colour Recognition %
                            
h = fspecial('Gaussian',[3 3],2); %Creates a Gaussian filter. 
blurred = imfilter(transformedImg,h); %Denoises the image using the Gaussian filter.
%figure,imshow(blurred)
se = strel('square',11); %Creates a structure element, that will be used later.
erosion = imerode(blurred,se); %Erosion performed,using the structure element to clear the blurred image.
%figure,imshow(erosion)
dilation = imdilate(erosion,se); %Dilation followed for last noise removal.
figure,imshow(dilation),title('Filtered Transformed Image')
C = makecform('srgb2lab'); %Converts the image into the Lab color space. 
transformedImgLab = applycform(dilation, C);
%figure,imshow(transformedImgLab)
Gray = im2gray(transformedImgLab ); %Converts the image to the grayscale.
%figure,imshow(binary)
edges = edge(Gray,'canny'); %Finds the most important edges, using the "canny" method.
%figure,imshow(edges)
se1 = strel('square',3); %Creates a new structure element
sharpened = imcomplement(imdilate(edges,se1)); %Dilation performed,sharpening the edges. Also imcomplement used, reversing the black and white colours.
%figure,imshow(sharpened)
                                 % Calculating stats %
                                       
                                       
 CC = bwconncomp(sharpened);
 squareMeasurements = regionprops(sharpened,'Area');
 Areas3 = [squareMeasurements.Area];
 %Areas3
                      
                               % Isolating the Squares %
 
isolatedSquares=zeros(CC.ImageSize);
for p=1:CC.NumObjects  %loop through each image
    if squareMeasurements(p).Area < 6500 && squareMeasurements(p).Area > 4000
        isolatedSquares(CC.PixelIdxList{p}) = 1; %set the image
    end

end

%figure,imshow(isolatedsquares,[0 1])
cleaned = imdilate(isolatedSquares,ones(7)); %Erosion is taking place to fill the last gaps.
%figure,imshow(cleaned), title('Isolated Squares')

       % Investigating the Centroid Coordinates of the isolated squares %
                           
                       
labeledSquareImg = bwlabel(cleaned);
squareMeasurements = regionprops(labeledSquareImg, 'Centroid','PixelIdxList');
allsquareCentroids = [squareMeasurements.Centroid];
%Centroids2=reshape(squareMeasurements,2,[]); 
%tr = transpose(Centroids2)

% Manually recording the Centroid Coordinates of each square in the images.
centroids =[ 113 113; 113 198 ; 113 283; 113 368; 198 113 ; 198 198 ; 198 283; 198 368; 283 113 ; 283 198; 283 283; 283 368; 368 113 ; 368 198; 368 283 ; 368 368];

                 
                      % Colours Translation into their L,a,b values %

white = [245,127,128]; 
blue = [81,145,61]; 
yellow = [231,114,214];
green = [227,74,204]; 
red = [133,190,192];


colors = cell(4);

for i = 1:16 %Loops through each square
    
    xmin = centroids(i); %x value of the centroid
    ymin = centroids(i,2); %y value of the centroid
    I2 = imcrop(transformedImgLab,[xmin ymin 15 15]); %Crops a small window around the centroid coordinates of each square with width,height = 15
    L_Channel = I2(:,:,1); %L values
    a_Channel = I2(:,:,2); %a values
    b_Channel = I2(:,:,3); %b values
    meanL = mean(L_Channel); %Calculates the mean value
    meanA = mean(a_Channel);
    meanB = mean(b_Channel);
    thisColor = [meanL(i),meanA(i), meanB(i)]; %Records the mean Lab values, found above in an array.
    distanceToWhite = sqrt((thisColor(1) - white(1)) .^ 2 + (thisColor(2) - white(2)) .^ 2 + (thisColor(3) - white(3)) .^ 2); %Calculates the Euclidean distance between the fixed colours and the colour that is being processed at a time.
    distanceToYellow = sqrt((thisColor(1) - yellow(1)) .^ 2 + (thisColor(2) - yellow(2)) .^ 2+ (thisColor(3) - yellow(3)) .^ 2);
    distanceToRed = sqrt((thisColor(1) - red(1)) .^ 2 + (thisColor(2) - red(2)) .^ 2 + (thisColor(3) - red(3)) .^ 2 );
    distanceToBlue = sqrt((thisColor(1) - blue(1)) .^ 2 + (thisColor(2) - blue(2)) .^ 2 + (thisColor(3) - blue(3)) .^ 2 );
    distanceToGreen = sqrt((thisColor(1) - green(1)) .^ 2 + (thisColor(2) - green(2)) .^ 2 + (thisColor(3) - green(3)) .^ 2);
    distances = [distanceToWhite,distanceToYellow,distanceToRed,distanceToBlue,distanceToGreen]; %Records all the distances to an array. 
    if distanceToWhite == min(distances) %Conditional statement that checks which fixed colour has the minimum distance to the colour that is being processed at that time. Then it appends the matching colour into the cell array and so on.
        colors(i) = cellstr('W'); 
    elseif distanceToYellow == min(distances)
        colors(i) = cellstr('Y');
    elseif distanceToRed == min(distances)
        colors(i) = cellstr('R');
    elseif distanceToBlue == min(distances)
        colors(i) = cellstr('B');
    elseif distanceToGreen== min(distances)
        colors(i) = cellstr('G');
    end
end

Results = reshape(colors,[],4);
result = Results
end
