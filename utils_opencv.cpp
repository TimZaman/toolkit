//utils_opencv.cpp

#include "utils_opencv.h"



bool util::detectClipping(cv::Mat matImage, int threshold_min, double percent_allowed_min, int threshold_max, double percent_allowed_max, std::string &strError){
	//cout << "detectClipping()" << endl;

	cv::Mat matGray;
	cvtColor(matImage, matGray, cv::COLOR_BGR2GRAY);
	if (matGray.depth()==CV_16U){
		//Convert to 8-bit
		matGray.convertTo(matGray, CV_8U, 1.0/257.0);
	}

	int numMin=0;
	int numMax=0;
	int numTotal = matGray.rows*matGray.cols;
	for (int i=0; i<numTotal; i++){
		int val = (int)matGray.at<unsigned char>(i);
		if (val <= threshold_min){
			numMin++;
		}
		if (val >= threshold_max){
			numMax++;
		}
	}

	double percent_min = 100*double(numMin)/double(numTotal);
	double percent_max = 100*double(numMax)/double(numTotal);

	//cout << " values < " << threshold_min << " = " << percent_min << "\%" << endl;
	//cout << " values > " << threshold_max << " = " << percent_max << "\%" << endl;

	if (percent_min > percent_allowed_min){
		strError = "Underexposure detected: " + std::to_string(percent_min) + "\% of values are below or equal to " + std::to_string(threshold_min) + "\n";
		return false;
	}
	if (percent_max > percent_allowed_max){
		strError = "Overexposure detected: " + std::to_string(percent_max) + "\% of values are above or equal to " + std::to_string(threshold_max) + "\n";
		return false;
	}
	return true;
}



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



std::vector<std::vector<int> > util::groupPoints(std::vector<cv::Point> vecPts, double mindist, int mingroupsize){
	//Groups points together through recursive chaining, using euclidean distance as metric
	//Any point can only belong to one group
	//mindist [px] minimal distance to any other group member
	std::vector<std::vector<int> > vecPtGroups;

	//std::cout << "vecPts.size()=" << vecPts.size() << std::endl;

	std::vector< std::vector<int> > vecCloseTo(vecPts.size());
	//For each point, we make its own array that includes all the indices to the points that it's close to.
	for (int i=0; i<vecPts.size();i++){
		vecCloseTo[i].push_back(i); //Itself
		//for (int j=i+1; j<vecPts.size();j++){
		for (int j=0; j<vecPts.size();j++){
			if (i==j) continue;
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

	std::vector<bool> isGrouped(vecCloseTo.size(), 0);

	//NEW!
	for (int i=0; i<vecCloseTo.size();i++){
		//std::cout << "#G" << i << " : ";

		if (isGrouped[i]) continue;

		//Copy the current array to the newly formed group
		std::vector<int> vecGroup(vecCloseTo[i]);

		//Make the initial dynamic expanding group
		for (int j=0; j<vecGroup.size(); j++){ //Vector size of vecGroup will get reevaluated each iteration
			//Evaluate current dependency's depencencies and add those to the end of the group, only if it's not a previous dependency.
			int depnow = vecGroup[j];
			//std::cout << "d(" << depnow << ")[";

			//If we havent used this dependencies' dependencies already, add those
			//Append current deps' dependencies.
			for (int k=0; k<vecCloseTo[depnow].size(); k++){ //Walk the deps
				//Check if it's already present
				int subdep = vecCloseTo[depnow][k];
				if ( std::find(vecGroup.begin(), vecGroup.end(), subdep) != vecGroup.end() ){
					//It's already in the group
				} else {
					//Not in group yet! add it.	
					vecGroup.push_back(subdep);	
					//std::cout << subdep << " ";
				}
			}
			///std::cout << "] ";
		}

		for (int j=0;j<vecGroup.size();j++){
			isGrouped[vecGroup[j]] = true;
		}

		vecPtGroups.push_back(vecGroup);
		//std::cout << std::endl;
	}

	//for (int i=0; i<vecPtGroups.size();i++){
	//	std::cout << "#GROUP " << i << " : ";
	//	for (int j=0; j<vecPtGroups[i].size(); j++){
	//		std::cout << vecPtGroups[i][j] << " ";
	//	}
	//	std::cout << std::endl;
	//}

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

	//Get it as close to zero as possible
	if (rRectNow.angle < -90.) {
		rRectNow.angle += 180.0;
	} else  if (rRectNow.angle > 90.) {
		rRectNow.angle -= 180.0;
	}

	if (rRectNow.angle < -45.) {
		rRectNow.angle += 90.0;
		std::swap(rRectNow.size.width, rRectNow.size.height);
	} else  if (rRectNow.angle > 45.) {
		rRectNow.angle -= 90.0;
		std::swap(rRectNow.size.width, rRectNow.size.height);
	}
	return rRectNow;
}


void util::expand(cv::Rect & rBounding, double pixels){
	//Expands a rectangle
	double hpx = pixels*0.5;
	rBounding = rBounding + cv::Size(pixels,pixels);
	rBounding = rBounding - cv::Point(hpx,hpx);
}


cv::Rect util::constrainRectInSize(cv::Rect rRect, cv::Size sImage){
	cv::Rect rImage(cv::Point(0,0), sImage);
	cv::Rect rIntersection =  rRect & rImage;
	return rIntersection;
}


cv::Mat util::crop(cv::Mat matImage, cv::RotatedRect rRect){

//	util::rectangle(matImage, rRect, cv::Scalar(255,0,0), 2);
//	circle(matImage, rRect.center, 1, cv::Scalar(255,0,0));

	cv::RotatedRect rRectNow = fixRotatedRect(rRect);

//	util::rectangle(matImage, rRect, cv::Scalar(0,0,255), 1);

//	imwrite("/Users/tzaman/Desktop/matImage.jpg", matImage);

	//It's a rotatedrect, so first get the bounding box

	cv::Rect rBounding = rRectNow.boundingRect();

	//Now fix the bounding rectangle
	if (rRectNow.size.width > rBounding.width){
		std::cout << "width of crop wider than width of bb" << std::endl;
		double d =rRectNow.size.width - rBounding.width; //Check out the difference
		expand(rBounding, d);
	}

	if (rRectNow.size.height > rBounding.height){
		std::cout << "height of crop wider than width of bb" << std::endl;
		double d =rRectNow.size.height - rBounding.height; //Check out the difference
		expand(rBounding, d);
	}

	std::cout << rBounding.tl() << std::endl;
	std::cout << rBounding.br() << std::endl;

	//Constrain it
	cv::Rect rBoundingInside = constrainRectInSize(rBounding, matImage.size());

	//Crop this off, and clone (because we rotate later)
	cv::Mat matBound = matImage(rBoundingInside).clone();
//	imwrite("/Users/tzaman/Desktop/matBound.jpg", matBound);

	//Not pad it if needed
	int pad_top   = rBounding.tl().y < 0 ? -rBounding.tl().y : 0;
	int pad_left  = rBounding.tl().x < 0 ? -rBounding.tl().x : 0;
	int pad_right = rBounding.br().y >= matImage.rows ? rBounding.br().y-matImage.rows-1 : 0;
	int pad_bot   = rBounding.br().x >= matImage.cols ? rBounding.br().x-matImage.rows-1 : 0;

	//std::cout << pad_top << " " << pad_left << " " << pad_right << " " << pad_bot << std::endl;

	if (pad_top > 0 || pad_left > 0 || pad_right > 0 || pad_bot > 0 ){
		copyMakeBorder(matBound, matBound, pad_top, pad_bot, pad_left, pad_right, cv::BORDER_CONSTANT, cv::Scalar(0,0,0));
//		imwrite("/Users/tzaman/Desktop/matBoundBorder.jpg", matBound);
	}

	//Now move the rotatedrect back by the part we have cropped off due to the bounding rect
	rRectNow.center = rRectNow.center - cv::Point2f(rBounding.x, rBounding.y); //add +0.5? or -0.5? or?

	//float angle = rRectNow.angle;

	//Rotate around center
	rotate(matBound, rRectNow.angle, matBound);
//	imwrite("/Users/tzaman/Desktop/matBoundRot.jpg", matBound);


	//Now we can crop, outward from the middle, with the size of the rotatedrect
	cv::Mat matCrop;
	getRectSubPix(matBound, rRectNow.size, rRectNow.center, matCrop);

//	imwrite("/Users/tzaman/Desktop/matCrop.jpg", matCrop);

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



