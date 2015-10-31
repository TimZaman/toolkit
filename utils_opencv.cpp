//utils_opencv.cpp

#include "utils_opencv.h"



std::string util::matToJpgString(cv::Mat matImage){
	//cout << "MatMeta::matToJpgString()" << endl;

	//First check if i actually have data assigned.
	if (matImage.data == NULL){
		std::cout << "No image data." << std::endl;
		return "";
	}

	std::vector<uchar> ubuf;
	std::vector<int> properties;
	imencode(".jpg", matImage, ubuf, properties);

	std::string str(ubuf.begin(), ubuf.end());

	return str;
}


cv::Rect util::constrainRectInSize(cv::Rect rCrop, cv::Size sImage){
	cv::Rect rImage(cv::Point(0,0), sImage);
	cv::Rect rIntersection =  rCrop & rImage;
	return rIntersection;
}


cv::RotatedRect util::scale(cv::RotatedRect rRect, double scale){
	return cv::RotatedRect(rRect.center*scale, cv::Size2f(rRect.size.width*scale, rRect.size.height*scale), rRect.angle);
}


cv::Rect util::retainCenterBlob(cv::Mat & matImage, int white_threshold){
	cv::Rect bounding_rect;
	//For 'Vrijstaand maken'
	cv::Mat matSmall;
	double smallfact = 0.1; //formerly 0.2
	cv::resize(matImage, matSmall, cv::Size(), smallfact, smallfact, cv::INTER_AREA);
	cv::cvtColor(matSmall, matSmall, cv::COLOR_BGR2GRAY);

	//imwrite("/Users/tzaman/Desktop/test.tif",matSmall);

	cv::Mat matThres;
	int thres_option = cv::THRESH_BINARY_INV; //INV for white
	//int white_threshold = this->autocrop_threshold; //just use the same autocrop toggle
	cv::threshold(matSmall, matThres, white_threshold, 255, thres_option);

	//imwrite("/Users/tzaman/Desktop/matThres.tif",matThres);

	//int largest_area=0;
	//int largest_contour_index=0;
	bool foundBlob=false;
	std::vector< std::vector<cv::Point> > contours; // Vector for storing contour
	std::vector< cv::Vec4i> hierarchy;
	cv::findContours( matThres, contours, hierarchy, CV_RETR_CCOMP, CV_CHAIN_APPROX_SIMPLE ); // Find the contours in the image
	for( int i = 0; i< contours.size(); i++ ) {// iterate through each contour. 
		/*double a=contourArea( contours[i],false);  //  Find the area of contour
		if(a>largest_area){
			largest_area=a;
			largest_contour_index=i;                //Store the index of largest contour
			bounding_rect=boundingRect(contours[i]); // Find the bounding rectangle for biggest contour
		}*/

		double pointInside = cv::pointPolygonTest(contours[i], cv::Point2f(matThres.cols*0.5, matThres.rows*0.5), false);
		if (pointInside>0){//Is inside
			bounding_rect=cv::boundingRect(contours[i]); // Find the bounding rectangle for biggest contour
			foundBlob=true;
			//break;
		}
	}
	//if (foundBlob){
	//	cout << "Found a blob at:" <<  bounding_rect << endl;
	//
	//} else {
	//	cout << WARNING << "Did NOT find blob. returning." << endl;
	//}

	return bounding_rect;
}


cv::Rect util::findBiggestBlob(cv::Mat & matImage){
	int largest_area=0;
	int largest_contour_index=0;
	cv::Rect bounding_rect(cv::Point(0,0), matImage.size());
	std::vector< std::vector<cv::Point> > contours; // Vector for storing contour
	std::vector< cv::Vec4i> hierarchy;
	cv::findContours( matImage, contours, hierarchy, CV_RETR_CCOMP, CV_CHAIN_APPROX_SIMPLE ); // Find the contours in the image
	for( int i = 0; i< contours.size(); i++ ) {// iterate through each contour. 
		double a=contourArea( contours[i],false);  //  Find the area of contour
		if(a>largest_area){
			largest_area=a;
			largest_contour_index=i;                //Store the index of largest contour
			bounding_rect=cv::boundingRect(contours[i]); // Find the bounding rectangle for biggest contour
		}
	}
	//drawContours( matImage, contours, largest_contour_index, Scalar(255), CV_FILLED, 8, hierarchy ); // Draw the largest contour using previously stored index.
	return bounding_rect;
}

cv::RotatedRect util::findBiggestBlobRot(cv::Mat & matImage){
	int largest_area=0;
	int largest_contour_index=0;
	cv::RotatedRect bounding_rotatedRect(cv::Point(0,0), matImage.size(),0);
	std::vector< std::vector< cv::Point> > contours; // Vector for storing contour
	std::vector< cv::Vec4i> hierarchy;
	cv::findContours( matImage, contours, hierarchy, CV_RETR_CCOMP, CV_CHAIN_APPROX_SIMPLE ); // Find the contours in the image
	for( int i = 0; i< contours.size(); i++ ) {// iterate through each contour. 
		double a=cv::contourArea( contours[i],false);  //  Find the area of contour
		if(a>largest_area){
			largest_area=a;
			largest_contour_index=i;                //Store the index of largest contour
			bounding_rotatedRect=cv::minAreaRect(contours[i]);
		}
	}
	//drawContours( matImage, contours, largest_contour_index, Scalar(255), CV_FILLED, 8, hierarchy ); // Draw the largest contour using previously stored index.
	return bounding_rotatedRect;
}



cv::Point2f util::ptMove(cv::Point2f pt, double dist, double angle_deg){
	double radians = angle_deg * M_PI / 180.0;
	pt.x = pt.x + dist * cos(radians);
	pt.y = pt.y + dist * sin(radians);	
	return pt;
}


cv::Point util::rect2cp(cv::Rect rRect){
	//Finds the center point of a rectangle
	return cv::Point(rRect.x+rRect.width*0.5, rRect.y+rRect.height*0.5);
}


std::vector<cv::Point> util::vecrotrect2vecpt(std::vector<cv::RotatedRect> vecRotRect){
	//Converts an array of RotatedRect's to an array of their center points.
	std::vector<cv::Point> vecPts(vecRotRect.size());
	for(int i=0; i<vecRotRect.size(); i++){
		vecPts[i] = cv::Point(vecRotRect[i].center.x, vecRotRect[i].center.y);
	}
	return vecPts;
}


void util::addRecursive(std::vector<int> & group, int myid, std::vector< std::vector<int> > & vecCloseTo, std::vector<int> & alreadyInGroup){
	//std::cout << "addRecursive(" << myid << "..)" << std::endl;
	//'myid' is the current index we are iterating through

	//Skip if i was already in group
	for (int i=0; i<alreadyInGroup.size(); i++){
		if (myid==alreadyInGroup[i]) {
			return;
		}
	}
	alreadyInGroup.push_back(myid); //Add myself to group

	//std::cout << myid << " will add [";
	//for (int i=0; i < vecCloseTo[myid].size(); i++){
	//	std::cout << vecCloseTo[myid][i] <<" ";
	//}
	//std::cout << "]" << std::endl;


	for (int i=0; i < vecCloseTo[myid].size(); i++){
		int depnow = vecCloseTo[myid][i];
		for (int i=0; i<alreadyInGroup.size(); i++){
			if (depnow==alreadyInGroup[i]) {
				continue;
			}
		}
		group.push_back(depnow);
		addRecursive(group, depnow, vecCloseTo, alreadyInGroup);
	}
}



std::vector<std::vector<int> > util::groupPoints(std::vector<cv::Point> vecPts, double mindist, int mingroupsize){
	//Groups points together through recursive chaining, using euclidean distance as metric
	//Any point can only belong to one group
	//mindist [px] minimal distance to any other group member
	std::vector<std::vector<int> > vecPtGroups;

	std::cout << "vecPts.size()=" << vecPts.size() << std::endl;

	std::vector< std::vector<int> > vecCloseTo(vecPts.size());
	//First make a list of stuff anything is close to withing the distance
	for (int i=0; i<vecPts.size();i++){
		//vecCloseTo[i].push_back(i); //Itself
		for (int j=i+1; j<vecPts.size();j++){
		//for (int j=0; j<vecPts.size();j++){
			double dist = util::pointDist(vecPts[i], vecPts[j]);
			if (dist<mindist){
				vecCloseTo[i].push_back(j);
			}
		}
	}

	//for (int i=0; i<vecCloseTo.size();i++){
	//	std::cout << "#" << i << " : ";
	//	for (int j=0; j<vecCloseTo[i].size(); j++){
	//		std::cout << vecCloseTo[i][j] << " ";
	//	}
	//	std::cout << std::endl;
	//}

	std::vector<int> alreadyInGroup;
	for (int i=0; i<vecCloseTo.size();i++){	
		for (int j=0; j<alreadyInGroup.size();j++){
			if (i==alreadyInGroup[j]) continue;
		}

		//std::cout << std::endl << " = New Group = " << std::endl;
		//Now we have to use an self-recursive 'Escher' statement to follow the rabbit down the hole
		std::vector<int> group;
		addRecursive(group, i, vecCloseTo, alreadyInGroup);

		//std::cout << "Group now:" << std::endl;
		//for (int j=0; j<group.size(); j++){
		//	std::cout << group[j] << " ";
		//}
		//std::cout << std::endl;


		if (group.size() > mingroupsize){
			//Remove duplicates
			std::sort( group.begin(), group.end());
			group.erase( unique( group.begin(), group.end() ), group.end() );
			vecPtGroups.push_back(group);
		}
	}

	std::cout << "vecPtGroups.size()=" << vecPtGroups.size() << std::endl;
	for (int i=0; i<vecPtGroups.size(); i++){
		std::cout << "vecPtGroups[" << i << "].size()=" << vecPtGroups[i].size() << std::endl;	
	}

	std::cout << "END groupPoints()" << std::endl;
	return vecPtGroups;
}





void util::rectangle(cv::Mat matImage, cv::RotatedRect rRect, cv::Scalar color, int thickness){
	//Draws a rectangle
	cv::Point2f rect_points[4]; 
	rRect.points( rect_points );
	for( int j = 0; j < 4; j++ ) {
		line( matImage, rect_points[j], rect_points[(j+1)%4], color, thickness, 8 );
	}
}



cv::Point2f util::rotate2d(const cv::Point2f& inPoint, const double& angRad)
{
    cv::Point2f outPoint;
    //CW rotation
    outPoint.x = std::cos(angRad)*inPoint.x - std::sin(angRad)*inPoint.y;
    outPoint.y = std::sin(angRad)*inPoint.x + std::cos(angRad)*inPoint.y;
    return outPoint;
}

cv::Point2f util::rotatePoint(const cv::Point2f& inPoint, const cv::Point2f& center, const double& angRad)
{
    return rotate2d(inPoint - center, angRad) + center;
}

double util::pts2angleDeg(cv::Point pt1, cv::Point pt2){
	double angledeg=atan2(double(pt2.y-pt1.y),double(pt2.x-pt1.x))*180.0/M_PI;
	return angledeg;
	//      -90
	// -180 \|/  -0
	//      -o-
	//  +0  /|\   180
	//      +90
}


cv::RotatedRect util::fixRotatedRect(cv::RotatedRect rRect){
	//Fixes a large angular rotation to a flip in width
	cv::RotatedRect rRectNow = rRect;
	if (rRectNow.angle < -45.) {
		rRectNow.angle += 90.0;
		std::swap(rRectNow.size.width, rRectNow.size.height);
	} else  if (rRectNow.angle > 45.) {
		rRectNow.angle -= 90.0;
		std::swap(rRectNow.size.width, rRectNow.size.height);
	}
	return rRectNow;
}

//OpenCV related
cv::Mat util::crop(cv::Mat matImage, cv::RotatedRect rRect){

	cv::RotatedRect rRectNow = fixRotatedRect(rRect);

	//It's a rotatedrect, so first get the bounding box

	cv::Rect rBounding = rRectNow.boundingRect();

	//Constrain it
	rBounding = constrainRectInSize(rBounding, matImage.size());

	//Crop this off, and clone (because we rotate later)

	cv::Mat matBound = matImage(rBounding).clone();

	//Now move the rotatedrect back by the part we have cropped off due to the bounding rect
	rRectNow.center = rRectNow.center - cv::Point2f(rBounding.x, rBounding.y); //add +0.5? or -0.5? or?

	float angle = rRectNow.angle;

	//Rotate around center
	rotate(matBound, angle, matBound);
	//imwrite("matBoundRot.png", matBound);

	cv::Size rect_size = rRectNow.size;
	
	//Account for rotation
	if (rRectNow.angle < -45.) {
		angle += 90.0;
		std::swap(rect_size.width, rect_size.height);
	} else  if (rRectNow.angle > 45.) {
		angle -= 90.0;
		std::swap(rect_size.width, rect_size.height);
	}


	//Now we can crop, outward from the middle, with the size of the rotatedrect
	cv::Mat matCrop;
	getRectSubPix(matBound, rect_size, rRectNow.center, matCrop);

	return matCrop;

}


void util::autoClipBrighten(cv::Mat &matImage, double percentile_lower, double percentile_upper){
	//First, auto brighten the crop. Find the brightness of the 5% percentile and 95%, and clip those to [0-255]

	//Check parameters
	if (percentile_lower <0 || percentile_upper <0){
		std::cout << "Error in autoClipBrighten: percentile is lower than 0." << std::endl;
		return;
	} else if (percentile_lower>1 || percentile_upper>1){
		std::cout << "Error in autoClipBrighten: percentile exceeds 1." << std::endl;
		return;
	} else if (percentile_lower >= percentile_upper){
		std::cout << "Error in autoClipBrighten: lower percentile lower than high one." << std::endl;
		return;
	}


	float range[] = { 0, 256 } ; //the upper boundary is exclusive
	const float* histRange = { range };
	bool uniform = true; bool accumulate = false;
	cv::Mat hist255;
	int histSize255 = 255;
	cv::calcHist( &matImage, 1, 0, cv::Mat(), hist255, 1, &histSize255, &histRange, uniform, accumulate ); //last arg=do accumulate
	//Normalize such that the area is 1 (not the peaks!)
	double histarea255=0;
	for (int i=0;i<255;i++){
	  histarea255 += hist255.at<float>(i);
	}
	//std::cout << "histarea255=" << histarea255 << std::endl;
	//Account for the histogram area, set cum to 1.
	float hist255cum=0;
	int binlow=-1; //lower percentile
	int binhigh=-1; //higher percentile
	for (int i=0;i<255;i++){
	  float valnow = hist255.at<float>(i) * (1.0/histarea255);
	  hist255cum+=valnow; //Add to cum
	  hist255.at<float>(i) = valnow;
	  if (hist255cum>percentile_lower && binlow==-1){
	    binlow=i;
	  }
	  if (hist255cum>percentile_upper && binhigh==-1){
	    binhigh=i;
	  }
	  //std::cout <<  "histcum=" << hist255cum << std::endl;
	}
	//std::cout << hist255;

	//std::cout << "binhigh=" << binhigh << " low=" << binlow << std::endl;

	//Now scale the image to the percentiles, 5%-95% to [0-255]
	//Caculate the offset and scaling values
	double ab_offset = -binlow;//auto brightness offset value
	//Then the scale value after that
	double ab_scale = 255.0/(binhigh-binlow);

	matImage = (matImage+ab_offset)*ab_scale;
	//matImage is passed by by reference
}


cv::Size util::getFitSize(cv::Size sizeIn, cv::Size sizeOut){
	//Gets the output size for 'fittin within'
	double dAspRat    = double(sizeIn.width)/double(sizeIn.height);  //Aspect ratio of the image itself
	double dAspRatBox = double(sizeOut.width)/double(sizeOut.height); //Aspect ratio of the box 
	cv::Size sizeFit;
	if (dAspRat>dAspRatBox){ //Match width
		sizeFit = cv::Size(sizeOut.width, floor(sizeOut.height*dAspRatBox/dAspRat));
	} else { //Match height
		sizeFit = cv::Size(floor(sizeOut.width*dAspRat/dAspRatBox), sizeOut.height);
	}
	return sizeFit;
}


double util::pointDist(cv::Point pt1, cv::Point pt2){
  return sqrt(pow((double)(pt1.x-pt2.x),2)+pow((double)(pt1.y-pt2.y),2));
}

double util::pointDist(cv::Point2f pt1, cv::Point2f pt2){
  return sqrt(pow((double)(pt1.x-pt2.x),2)+pow((double)(pt1.y-pt2.y),2));
}

void util::rot90(cv::Mat &matImage, int rotflag){
	//1=CW, 2=CCW, 3=180
	if (rotflag == 1){
		transpose(matImage, matImage);	
		flip(matImage, matImage,1);	//transpose+flip(1)=CW
	} else if (rotflag == 2) {
		transpose(matImage, matImage);	
		flip(matImage, matImage,0);	//transpose+flip(0)=CCW		
	} else if (rotflag ==3){
		flip(matImage, matImage,-1);	//flip(-1)=180			
	} else if (rotflag != 0){ //if not 0,1,2,3:
		std::cout  << "Unknown rotation flag(" << rotflag << ")" << std::endl;
	}
}

void util::rotate(cv::Mat& src, double angle, cv::Mat& dst){
    //std::cout << RANDCOL << "R O T A T I N G" << std::endlr;
    //int len = std::max(src.cols, src.rows);
    cv::Point2f ptCp(src.cols*0.5, src.rows*0.5);
    //cv::Point2f pt(len/2., len/2.);
    cv::Mat M = cv::getRotationMatrix2D(ptCp, angle, 1.0);
    cv::warpAffine(src, dst, M, src.size(), cv::INTER_CUBIC); //Nearest is too rough, 
}



