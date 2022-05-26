
//Author: Annalisa Bellandi
//Faulkner lab, John Innes Centre

///This macro allows for:
///creation of 10 rays perpendicular to a defined segment
///scanning fluorescence intensity along these rays over time

/// Main steps:
/// The users selects two points along the axis of interest
/// 10 rays are created perpendicular to the axis
/// fluorescence intensity along the rays is saved for each time step


//---------------------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------------

//BEFORE STRATING

//open iamge
//note stable frame (sf)
//note frame of when veins are visible (vv)
//run stackreg 
//proceed with script
    

//--------------------------------------------------------------------------------------------------------------------------

dir = getInfo("image.directory");// get the image directory
print(dir)
imageTitle=getTitle; //get the image name
					
//get pixel size
getPixelSize(unit, pixelWidth, pixelHeight);
print("Current image pixel width = " + pixelWidth + " " + unit +".");
pxw = pixelWidth;
print(pxw);
					
//check planeTimings.txt https://docs.openmicroscopy.org/bio-formats/5.7.3/users/imagej/
run("Bio-Formats Macro Extensions"); //to get metadata info 
selectWindow(imageTitle); //have to return to the original because I didn't save duplicate in folder
id = getInfo("image.directory") + getInfo("image.filename");
Ext.setId(id);
Ext.getImageCount(imageCount);
deltaT = newArray(nSlices); //create array with vector of time points
title_without_extension = substring(imageTitle, 0, lengthOf(imageTitle)-4);
					
close(); //close orignal
					
run("Bio-Formats Importer", "open=[" + dir + title_without_extension + "_reg.tif" + "] autoscale color_mode=Default open_files rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");

run("Duplicate...", "duplicate"); //duplicate the entire stack for measurements
measuredImageTitle="dup"+imageTitle; //add the starting "dup" to the original image title
rename(measuredImageTitle);//rename it
selectWindow(measuredImageTitle); //select the duplicate image to measure

run("Clear Results");

waitForUser("make point at the start of the future line");
getSelectionCoordinates(x, y);
 xa = x[0];
 ya = y[0];

waitForUser("make intermediate point, where you want your line to pass");
getSelectionCoordinates(x, y);
 xb = x[1];
 yb = y[1];

H = yb - ya
L = xb - xa

makeLine(xa,ya,xb,yb);

alfa = atan(H/L)//*180/PI the formulas take radiants so no need to convert in degrees
print(alfa)

//set some parameters
K = 200/pxw //calculates how long is the line along the vein in pixels depending on how long you want it in microns
ntransect = 10 // how many transects you want to have
stepsize = K/ntransect // what is the distacne between each transect
widthlines = 10/pxw //calculates the width of the line is pixels depending on how wide you want your line in micron
S = 500/pxw //S is the leght of each transect

xc = xa+K*cos(alfa)
yc = ya+K*sin(alfa)

//this is the line along the vein
makeLine(xa,ya,xc,yc);
Overlay.addSelection("yellow", 1);

beta = alfa + (PI/2)
print(beta)
gamma = PI - beta

//Loop0: repeats the loop1 and 2 for each time point (slice) - the loops starts like this: for(k=0; k<=nSlices; k++){ , insert slice in place of 0
for(k=1; k<=nSlices; k++){
	setSlice(k);
	no=k-1; // arrays start at zero, so I subtract 1.
	Ext.getPlaneTimingDeltaT(deltaT[no], no); //....
	t=deltaT[no];
	print(t);
    
//loop1: for each step along the chosen vein section, creates a perfendicular line and gets the profile
  for(i=0; i<=K; i+=stepsize){
	 x1 = xa+i*cos(alfa);
	 y1 = ya+i*sin(alfa);
	 makeLine(x1,y1, x1-(S*cos(gamma)), y1+(S*sin(gamma)),widthlines);
	 Overlay.addSelection("green", widthlines);
     profile = getProfile();
  
//Loop2: for each  line perpendicular to the vein, saves the values of the profile and the position of the values along the lenght in a table
     for(j=0; j<profile.length; j++){
    	 setResult(i*pxw, j+S*(k-1), profile[j]);
    	 setResult("frame", j+S*(k-1), k);
    	 setResult("t", j+S*(k-1), t); 
    	 //once you got the time, put t as a value in place of k: setResult("t", j+50*(k-1), t);
    	 setResult("d", j+S*(k-1), j*pxw); // this should give the space in um in my d column instead fo the pxl n
       
        }
    }
}

title_without_extension = substring(imageTitle, 0, lengthOf(imageTitle)-4);
saveAs("Results", dir + title_without_extension + ".csv");

waitForUser("Happy with the scan?");
close();


