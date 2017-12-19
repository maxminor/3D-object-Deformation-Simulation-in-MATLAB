myPde = createpde(1);
importGeometry(myPde,'./obj/teapot.obj');
figure
pdegplot(myPde,'FaceLabels','on');