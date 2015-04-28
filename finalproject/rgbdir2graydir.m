%convert a directory of color images to gray images
imdir = dir('08/image_2');
for i=0:length(imdir)-2
    I = imread(['08/image_2/' num2str(i,'%06d') '.png']);
    I = rgb2gray(I);
    imwrite(I, ['08/mono_gray/' num2str(i,'%06d') '.png']);
end