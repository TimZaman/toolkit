//utils_opencv.cpp

#include "utils_opencv.h"


//OpenCV related
cv::Mat util::crop(cv::Mat matImage, cv::RotatedRect rRect){


	//It's a rotatedrect, so first get the bounding box

	cv::Rect rBounding = rRect.boundingRect();

	//Constrain it
	rBounding = constrainRectInSize(rBounding, matImage.size());

	//Crop this off, and clone (because we rotate later)

	cv::Mat matBound = matImage(rBounding).clone();

	//Now move the rotatedrect back by the part we have cropped off due to the bounding rect
	//cv::RotatedRect rRectSmall = rRect - Point2f(2,2);//rBounding.tl();
	
	rRect.center = rRect.center - cv::Point2f(rBounding.x, rBounding.y); //add +0.5? or -0.5? or?

	float angle = rRect.angle;

	//Rotate around center
	rotate(matBound, angle, matBound);
	//imwrite("matBoundRot.png", matBound);

	cv::Size rect_size = rRect.size;
	
	//Account for rotation
	if (rRect.angle < -45.) {
		angle += 90.0;
		swap(rect_size.width, rect_size.height);
	} else  if (rRect.angle > 45.) {
		angle -= 90.0;
		swap(rect_size.width, rect_size.height);
	}


	//Now we can crop, outward from the middle, with the size of the rotatedrect
	cv::Mat matCrop;
	getRectSubPix(matBound, rect_size, rRect.center, matCrop);

	return matCrop;

}


void util::autoClipBrighten(cv::Mat &matImage, double percentile_lower, double percentile_upper){
	//First, auto brighten the crop. Find the brightness of the 5% percentile and 95%, and clip those to [0-255]

	//Check parameters
	if (percentile_lower <0 || percentile_upper <0){
		cout << "Error in autoClipBrighten: percentile is lower than 0." << endl;
		return;
	} else if (percentile_lower>1 || percentile_upper>1){
		cout << "Error in autoClipBrighten: percentile exceeds 1." << endl;
		return;
	} else if (percentile_lower >= percentile_upper){
		cout << "Error in autoClipBrighten: lower percentile lower than high one." << endl;
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
	//cout << "histarea255=" << histarea255 << endl;
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
	  //cout <<  "histcum=" << hist255cum << endl;
	}
	//cout << hist255;

	//cout << "binhigh=" << binhigh << " low=" << binlow << endl;

	//Now scale the image to the percentiles, 5%-95% to [0-255]
	//Caculate the offset and scaling values
	double ab_offset = -binlow;//auto brightness offset value
	//Then the scale value after that
	double ab_scale = 255.0/(binhigh-binlow);

	matImage = (matImage+ab_offset)*ab_scale;
	//matImage is passed by by reference
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
		cout  << "Unknown rotation flag(" << rotflag << ")" << endl;
	}
}

void util::rotate(cv::Mat& src, double angle, cv::Mat& dst){
    //cout << RANDCOL << "R O T A T I N G" << endlr;
    //int len = std::max(src.cols, src.rows);
    cv::Point2f ptCp(src.cols*0.5, src.rows*0.5);
    //cv::Point2f pt(len/2., len/2.);
    cv::Mat M = cv::getRotationMatrix2D(ptCp, angle, 1.0);
    cv::warpAffine(src, dst, M, src.size(), cv::INTER_CUBIC); //Nearest is too rough, 
}

cv::Rect util::constrainRectInSize(cv::Rect rCrop, cv::Size sImage){
  cv::Rect rImage(cv::Point(0,0), sImage);
  cv::Rect rIntersection =  rCrop & rImage;
  return rIntersection;
}


