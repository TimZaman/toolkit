
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

    double percent_min = 100.0 * double(numMin) / double(numTotal);
    double percent_max = 100.0 * double(numMax) / double(numTotal);

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

double util::rectangleSimilarity(cv::Rect rect_lhs, cv::Rect rect_rhs) {
    cv::Rect rect_overlap = rect_lhs & rect_rhs;
    double overlap_area = static_cast<double>(rect_overlap.area());
    double overlap_lhs = overlap_area / static_cast<double>(rect_lhs.area());
    double overlap_rhs = overlap_area / static_cast<double>(rect_rhs.area());
    return std::min(overlap_lhs, overlap_rhs);
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
    //    cout << "Found a blob at:" <<  bounding_rect << endl;
    //
    //} else {
    //    cout << WARNING << "Did NOT find blob. returning." << endl;
    //}

    return bounding_rect;
}

cv::RotatedRect util::findRectAroundLargeBlobs(cv::Mat & mat_image, double minimum_blobfactor) {
    int largest_area = 0;
    int largest_contour_index = 0;

    int minimum_blob_area = mat_image.cols * mat_image.rows * minimum_blobfactor;

    std::vector < cv::RotatedRect > large_rotrects; 
    std::vector<cv::Point> point_on_large_contour;
    std::vector< std::vector<cv::Point> > contours; // Vector for storing contour
    std::vector< cv::Vec4i> hierarchy;

    cv::findContours(mat_image.clone(), contours, hierarchy, CV_RETR_CCOMP, CV_CHAIN_APPROX_SIMPLE); // Find the contours in the image
    for( int i = 0; i < contours.size(); i++ ) {
        double a = contourArea(contours[i], false);  //  Find the area of contour
        if (a > minimum_blob_area) {
            large_rotrects.push_back(cv::minAreaRect(contours[i]));
            point_on_large_contour.push_back(contours[i][0]);
        }
    }

    // Now draw lines between all rotrects, this connects all the largest blobs.
    for ( int i = 0; i < static_cast<int>(large_rotrects.size())-1; i++ ) {
        cv::line(mat_image, large_rotrects[i].center, large_rotrects[i+1].center, cv::Scalar(255), 2);
    }
    // Make sure the center line is connected with a point on the contour. This is allowed because i assume convexity
    for ( int i = 0; i < static_cast<int>(large_rotrects.size()); i++ ) {
        cv::line(mat_image, point_on_large_contour[i], large_rotrects[i].center, cv::Scalar(255), 2);
    }
    //cv::imwrite("/Users/tzaman/Desktop/mat_connected_blobrots.png", mat_image);

    // Finally we do a biggest blob (all big blobs are not connected! what a great dirty hack!)
    cv::RotatedRect bounding_rect = util::findBiggestBlobRot(mat_image);

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
            largest_contour_index = i;                //Store the index of largest contour
            bounding_rect = cv::boundingRect(contours[i]); // Find the bounding rectangle for biggest contour
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

cv::Point util::rect2cp(cv::Rect rotrect) {
    //Finds the center point of a rectangle
    return cv::Point(rotrect.x + rotrect.width * 0.5, rotrect.y + rotrect.height * 0.5);
}

std::vector<cv::Point> util::vecrotrect2vecpt(std::vector< cv::RotatedRect > vecRotRect) {
    //Converts an array of RotatedRect's to an array of their center points.
    std::vector<cv::Point> vecPts(vecRotRect.size());
    for(int i = 0; i < vecRotRect.size(); i++){
        vecPts[i] = cv::Point(vecRotRect[i].center.x, vecRotRect[i].center.y);
    }
    return vecPts;
}

std::vector<std::vector<int> > util::groupPoints(std::vector<cv::Point> vecPts, double mindist, int mingroupsize) {
    //Groups points together through recursive chaining, using euclidean distance as metric
    //Any point can only belong to one group
    //mindist [px] minimal distance to any other group member
    std::vector<std::vector<int> > vecPtGroups;

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

    std::vector<bool> isGrouped(vecCloseTo.size(), 0);

    for (int i = 0; i < vecCloseTo.size(); i++) {

        if (isGrouped[i]) continue;

        //Copy the current array to the newly formed group
        std::vector<int> vecGroup(vecCloseTo[i]);

        //Make the initial dynamic expanding group
        for (int j=0; j<vecGroup.size(); j++){ //Vector size of vecGroup will get reevaluated each iteration
            //Evaluate current dependency's depencencies and add those to the end of the group, only if it's not a previous dependency.
            int depnow = vecGroup[j];

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
                }
            }
        }

        for (int j=0;j<vecGroup.size();j++){
            isGrouped[vecGroup[j]] = true;
        }

        vecPtGroups.push_back(vecGroup);
    }

    return vecPtGroups;
}

void util::rectangle(cv::Mat image, cv::RotatedRect rotrect, cv::Scalar color, int thickness) {
    //Draws a rectangle
    cv::Point2f rect_points[4]; 
    rotrect.points( rect_points );
    for( int j = 0; j < 4; j++ ) {
        line( image, rect_points[j], rect_points[(j + 1) % 4], color, thickness, 8 );
    }
}

cv::Point2f util::rotate2d(const cv::Point2f& pt_in, const double& angle_rad) {
    cv::Point2f pt_out;
    //CW rotation
    pt_out.x = std::cos(angle_rad) * pt_in.x - std::sin(angle_rad) * pt_in.y;
    pt_out.y = std::sin(angle_rad) * pt_in.x + std::cos(angle_rad) * pt_in.y;
    return pt_out;
}

cv::Point2f util::rotatePoint(const cv::Point2f& pt_in, const cv::Point2f& center, const double& angle_rad) {
    return rotate2d(pt_in - center, angle_rad) + center;
}

double util::pts2angleRad(cv::Point pt1, cv::Point pt2) {
    return atan2(static_cast<double>(pt2.y - pt1.y), static_cast<double>(pt2.x - pt1.x));
    //      -?
    // -? \|/  -0
    //      -o-
    //  +0  /|\   ?
    //      ?
}

double util::pts2angleDeg(cv::Point pt1, cv::Point pt2) {
    return atan2(static_cast<double>(pt2.y - pt1.y), static_cast<double>(pt2.x - pt1.x)) * 180.0 / M_PI;
    //      -90
    // -180 \|/  -0
    //      -o-
    //  +0  /|\  +180
    //      +90
}

cv::RotatedRect util::fixRotatedRect(cv::RotatedRect rotrect) {
    //Fixes a large angular rotation to a flip in width
    cv::RotatedRect rotrect_fixed = rotrect;

    //Get it as close to zero as possible
    if (rotrect_fixed.angle < -90.) {
        rotrect_fixed.angle += 180.0;
    } else  if (rotrect_fixed.angle > 90.) {
        rotrect_fixed.angle -= 180.0;
    }

    if (rotrect_fixed.angle < -45.) {
        rotrect_fixed.angle += 90.0;
        std::swap(rotrect_fixed.size.width, rotrect_fixed.size.height);
    } else  if (rotrect_fixed.angle > 45.) {
        rotrect_fixed.angle -= 90.0;
        std::swap(rotrect_fixed.size.width, rotrect_fixed.size.height);
    }
    return rotrect_fixed;
}


void util::expand(cv::Rect & rect_bounding, double pixels) {
    //Expands a rectangle
    double hpx = pixels * 0.5;
    rect_bounding = rect_bounding + cv::Size(pixels, pixels) - cv::Point(hpx, hpx);
}


cv::Rect util::constrainRectInSize(cv::Rect rotrect, cv::Size size_image){
    cv::Rect rImage(cv::Point(0,0), size_image);
    cv::Rect rIntersection =  rotrect & rImage;
    return rIntersection;
}


cv::Mat util::crop(cv::Mat image, cv::RotatedRect rRect) {

    cv::RotatedRect rRectNow = fixRotatedRect(rRect);

    //It's a rotatedrect, so first get the bounding box

    cv::Rect rBounding = rRectNow.boundingRect();

    //Now fix the bounding rectangle
    if (rRectNow.size.width > rBounding.width) {
        std::cout << "width of crop wider than width of bb" << std::endl;
        double d =rRectNow.size.width - rBounding.width; //Check out the difference
        expand(rBounding, d);
    }

    if (rRectNow.size.height > rBounding.height) {
        std::cout << "height of crop wider than width of bb" << std::endl;
        double d =rRectNow.size.height - rBounding.height; //Check out the difference
        expand(rBounding, d);
    }

    //Constrain it
    cv::Rect rBoundingInside = constrainRectInSize(rBounding, image.size());

    //Crop this off, and clone (because we rotate later)
    cv::Mat matBound = image(rBoundingInside).clone();

    //Not pad it if needed
    int pad_top   = rBounding.tl().y < 0 ? - rBounding.tl().y : 0;
    int pad_left  = rBounding.tl().x < 0 ? - rBounding.tl().x : 0;
    int pad_right = rBounding.br().y >= image.rows ? rBounding.br().y - image.rows: 0;
    int pad_bot   = rBounding.br().x >= image.cols ? rBounding.br().x - image.cols: 0;

    //std::cout << " image =" << image.cols << "x" << image.rows << std::endl;
    //std::cout << " crop angle=" << rRect.angle << " center=" << rRect.center << " size=" << rRect.size << std::endl;
    //std::cout << " pad_top=" << pad_left << " pad_left=" << pad_left << " pad_right=" << pad_right << " pad_bot=" << pad_bot << std::endl;

    if (pad_top > 0 || pad_left > 0 || pad_right > 0 || pad_bot > 0 ) {
        copyMakeBorder(matBound, matBound, pad_top, pad_bot, pad_left, pad_right, cv::BORDER_CONSTANT, cv::Scalar(0,0,0));
    }

    //Now move the rotatedrect back by the part we have cropped off due to the bounding rect
    rRectNow.center = rRectNow.center - cv::Point2f(rBounding.x, rBounding.y); //add +0.5? or -0.5? or?

    //Rotate around center
    rotate(matBound, rRectNow.angle, matBound);

    //Now we can crop, outward from the middle, with the size of the rotatedrect
    cv::Mat matCrop;
    getRectSubPix(matBound, rRectNow.size, rRectNow.center, matCrop);

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
        flip(matImage, matImage,1);    //transpose+flip(1)=CW
    } else if (rotflag == 2) {
        transpose(matImage, matImage);    
        flip(matImage, matImage,0);    //transpose+flip(0)=CCW        
    } else if (rotflag ==3){
        flip(matImage, matImage,-1);    //flip(-1)=180            
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




static void util::rotatingCalipers( const cv::Point2f* points, int n, float* out ){
    // Function from OpenCV source code, changed/adapted from 'rotcalipers.cpp'
    //  cartesian coordinates are used
    float mindim = FLT_MAX;
    char buffer[32] = {};
    int i, k;
    cv::AutoBuffer<float> abuf(n*3);
    float* inv_vect_length = abuf;
    cv::Point2f* vect = (cv::Point2f*)(inv_vect_length + n);
    int left = 0, bottom = 0, right = 0, top = 0;
    int seq[4] = { -1, -1, -1, -1 };

    // Rotating calipers sides will always have coordinates:
    //  (a,b) (-b,a) (-a,-b) (b, -a)
    // 
    // This is a first base bector (a,b) initialized by (1,0)
    float orientation = 0;
    float base_a;
    float base_b = 0;

    float left_x, right_x, top_y, bottom_y;
    cv::Point2f pt0 = points[0];

    left_x = right_x = pt0.x;
    top_y = bottom_y = pt0.y;

    for (i = 0; i < n; i++) {
        double dx, dy;

        if( pt0.x < left_x ){
            left_x = pt0.x, left = i;
        }

        if( pt0.x > right_x ){
            right_x = pt0.x, right = i;
        }

        if( pt0.y > top_y ){
            top_y = pt0.y, top = i;
        }

        if( pt0.y < bottom_y ){
            bottom_y = pt0.y, bottom = i;
        }

        cv::Point2f pt = points[(i+1) & (i+1 < n ? -1 : 0)];

        dx = pt.x - pt0.x;
        dy = pt.y - pt0.y;

        vect[i].x = (float)dx;
        vect[i].y = (float)dy;
        inv_vect_length[i] = (float)(1./std::sqrt(dx*dx + dy*dy));

        pt0 = pt;
    }

    // Find convex hull orientation
    {
        double ax = vect[n-1].x;
        double ay = vect[n-1].y;

        for (i = 0; i < n; i++) {
            double bx = vect[i].x;
            double by = vect[i].y;

            double convexity = ax * by - ay * bx;

            if( convexity != 0 )
            {
                orientation = (convexity > 0) ? 1.f : (-1.f);
                break;
            }
            ax = bx;
            ay = by;
        }
        CV_Assert( orientation != 0 );
    }
    base_a = orientation;

    //****************************************************************************************/
    // Init calipers position
    seq[0] = bottom;
    seq[1] = right;
    seq[2] = top;
    seq[3] = left;
    //****************************************************************************************
    // Main loop - evaluate angles and rotate calipers

    // All of the edges will be checked while rotating calipers by 90 degrees
    for (k = 0; k < n; k++) {
        // compute cosine of angle between calipers side and polygon edge
        // dp - dot product
        float dp[4] = {
            +base_a * vect[seq[0]].x + base_b * vect[seq[0]].y,
            -base_b * vect[seq[1]].x + base_a * vect[seq[1]].y,
            -base_a * vect[seq[2]].x - base_b * vect[seq[2]].y,
            +base_b * vect[seq[3]].x - base_a * vect[seq[3]].y,
        };

        float maxcos = dp[0] * inv_vect_length[seq[0]];

        // Number of calipers edges, that has minimal angle with edge
        int main_element = 0;

        // Choose minimal angle
        for (i = 1; i < 4; ++i) {
            float cosalpha = dp[i] * inv_vect_length[seq[i]];
            if (cosalpha > maxcos) {
                main_element = i;
                maxcos = cosalpha;
            }
        }

        // Rotate Calipers
        {
            //get next base
            int pindex = seq[main_element];
            float lead_x = vect[pindex].x*inv_vect_length[pindex];
            float lead_y = vect[pindex].y*inv_vect_length[pindex];
            switch( main_element ) {
                case 0:
                    base_a = lead_x;
                    base_b = lead_y;
                    break;
                case 1:
                    base_a = lead_y;
                    base_b = -lead_x;
                    break;
                case 2:
                    base_a = -lead_x;
                    base_b = -lead_y;
                    break;
                case 3:
                    base_a = -lead_y;
                    base_b = lead_x;
                    break;
            default:
                CV_Error(CV_StsError, "main_element should be 0, 1, 2 or 3");
            }
        }
        // Change base point of main edge
        seq[main_element] += 1;
        seq[main_element] = (seq[main_element] == n) ? 0 : seq[main_element];

        // Find area of rectangle
        {
            float height;
            float longest_side;

            // find vector left-right
            float dx = points[seq[1]].x - points[seq[3]].x;
            float dy = points[seq[1]].y - points[seq[3]].y;

            // dotproduct */
            float width = dx * base_a + dy * base_b;

            // find vector left-right
            dx = points[seq[2]].x - points[seq[0]].x;
            dy = points[seq[2]].y - points[seq[0]].y;

            // dotproduct
            height = -dx * base_b + dy * base_a;

            longest_side = std::max(std::abs(width), std::abs(height));

            //Retain the one with the smallest long side
            if (longest_side <= mindim) {
                float *buf = (float *) buffer;

                mindim = longest_side;
                // leftist point
                ((int *) buf)[0] = seq[3];
                buf[1] = base_a;
                buf[2] = width;
                buf[3] = base_b;
                buf[4] = height;
                // bottom point
                ((int *) buf)[5] = seq[0];
            }
        }
    }

    {
        float *buf = (float *) buffer;

        float A1 = buf[1];
        float B1 = buf[3];

        float A2 = -buf[3];
        float B2 = buf[1];

        float C1 = A1 * points[((int *) buf)[0]].x + points[((int *) buf)[0]].y * B1;
        float C2 = A2 * points[((int *) buf)[5]].x + points[((int *) buf)[5]].y * B2;

        float idet = 1.f / (A1 * B2 - A2 * B1);

        float px = (C1 * B2 - C2 * B1) * idet;
        float py = (A1 * C2 - A2 * C1) * idet;

        out[0] = px;
        out[1] = py;

        out[2] = A1 * buf[2];
        out[3] = B1 * buf[2];

        out[4] = A2 * buf[4];
        out[5] = B2 * buf[4];
    }
}

cv::RotatedRect util::minAreaSquare( cv::InputArray _points ){
    //Adapted from OpenCV source code.
    //It is changed such that it gives the minimal enclosing *square* as opposed to rectangle.
    //It's very possible it will not actually yield a square shape, but you can expand
    //all sides to the widest one, this will be the enclosing square itself.
    cv::Mat hull;
    cv::Point2f out[3];
    cv::RotatedRect box;

    cv::convexHull(_points, hull, true, true);

    if( hull.depth() != CV_32F ) {
        cv::Mat temp;
        hull.convertTo(temp, CV_32F);
        hull = temp;
    }

    int n = hull.checkVector(2);
    const cv::Point2f* hpoints = hull.ptr<cv::Point2f>();

    if( n > 2 ) {
        rotatingCalipers( hpoints, n, (float*)out );
        box.center.x = out[0].x + (out[1].x + out[2].x)*0.5f;
        box.center.y = out[0].y + (out[1].y + out[2].y)*0.5f;
        box.size.width = (float)std::sqrt((double)out[1].x*out[1].x + (double)out[1].y*out[1].y);
        box.size.height = (float)std::sqrt((double)out[2].x*out[2].x + (double)out[2].y*out[2].y);
        box.angle = (float)atan2( (double)out[1].y, (double)out[1].x );
    } else if( n == 2 ) {
        box.center.x = (hpoints[0].x + hpoints[1].x)*0.5f;
        box.center.y = (hpoints[0].y + hpoints[1].y)*0.5f;
        double dx = hpoints[1].x - hpoints[0].x;
        double dy = hpoints[1].y - hpoints[0].y;
        box.size.width = (float)std::sqrt(dx*dx + dy*dy);
        box.size.height = 0;
        box.angle = (float)atan2( dy, dx );
    } else {
        if( n == 1 ) {
            box.center = hpoints[0];
        }
    }

    box.angle = (float)(box.angle*180/CV_PI);
    return box;
}

//#ifdef UTILS_GENERAL_TIM_H
cv::Mat util::correctGamma(cv::Mat &img, cv::Vec3d gamma , double contrast , bool forSaving){
    return correctGamma(img, gamma[2], gamma[1], gamma[0], contrast , forSaving);
}

cv::Mat util::correctGamma(cv::Mat &img, double gammaR,  double gammaG ,  double gammaB , double contrast , bool forSaving) {
    //cout << "inside correctGamma.." << endl;
    int depth;
    int maxbit;

    cv::Mat returnImg;

    if (forSaving == true){ //Leave bit-depth as-is for saving
        //cout << "image is for left as is.." << endl;
        img.copyTo(returnImg);
    } else {
        if (img.depth() == CV_8U){
            //cout << "image is for left as is (8bit).." << endl;
            img.copyTo(returnImg);
        } else if (img.depth() == CV_16U){
            //cout << "image is converted to 8 bit.." << endl;
            img.convertTo(returnImg, CV_8U, 1.0/257.0);
        }
    }


    if ( returnImg.depth() == CV_8U){
        //cout << "depth=CV_8U" << endl;
        depth = CV_8U;
        maxbit=255;
    } else if ( returnImg.depth() == CV_16U) {
        //cout << "depth=CV_16U" << endl;
        depth = CV_16U;
        maxbit = 65535;
    } else {
        std::cout << "WARNING! UNKNOWN BITDEPTH!" << std::endl;
        return img;
    }

    int N_SEG = maxbit+1;
    //cout << "N_SEG=" << N_SEG << endl;

    //Compute the difference between the scales it should be contained within
    //double dif = absMaxNow-absMinNow; //240-13=227



    //cout << "Correcting gamma.." << endl;
    //vector<Mat> imgSplit(returnImg.channels());
    std::vector<cv::Mat> splitResult(returnImg.channels());

    //cout << "splitting image.." << endl;
    cv::split(returnImg, splitResult);

    //double gammaArr[3] = {gammaR, gammaG, gammaB};
    double gammaArr[3] = {gammaB, gammaG, gammaR}; //OpenCV BGR style
    //Sanitize values  between [-1:1]
    contrast = (contrast < -0.5)? -0.5 : contrast;
    contrast = (contrast >  0.5)?  0.5 : contrast;

    //cout << "Starting loop.." << endl;

    //tic();
    //#pragma omp parallel for
    for (int c=0; c<returnImg.channels(); c++){

        double gamma = gammaArr[c];

        //Sanitize values  between [-1:1]
        gamma = (gamma < -0.5)? -0.5 : gamma;
        gamma = (gamma >  0.5)?  0.5 : gamma;


        //cout << "gamma=" << gamma << " contrast=" << contrast << endl;
        //imgSplit[c].copyTo(splitResult[c]);

        //double gamma = 0.3;
        //double contrast = 0.5;
        //cout << "making bezier.." << endl;
        //tic();

        int lutX[N_SEG];
        int lutY[N_SEG];    
        //std::cout << "making bezier.." << std::endl;
        util::makeBezier(gamma, contrast, N_SEG, lutX, lutY);
        //toc();
        //std::cout << "bezier done." << std::endl;

        //for (int i=0;i<N_SEG;i=i+1){
        //    cout << lutX[i] << "," << lutY[i] << endl;
        //}

        
        //std::cout << "applying LUT.." << std::endl;
        //tic();
        /*
        for(int x=0; x<returnImg.cols; x++){
            for(int y=0; y<returnImg.rows; y++){
                //8bit=Vec3b; 16bit=Vec3w
                //if (depth==CV_16U){
                //    returnImg.at<Vec3w>(y,x)[c] = lutY[ returnImg.at<Vec3w>(y,x)[c] ];
                //} else {
                //    returnImg.at<Vec3b>(y,x)[c] = lutY[ returnImg.at<Vec3b>(y,x)[c] 
                //}
            }
        } 
        */
        //toc();


        int dim(256);
        cv::Mat lut(1, &dim, CV_8U);
        for (int i=0; i<256; i++){
            //lut.at<uchar>(i)= 255-i;
            lut.at<uchar>(i)= lutY[i];
        }
        cv::LUT(splitResult[c], lut, splitResult[c]);

        //std::cout << "Done LUT application" << std::endl;
    }
    //toc();
    
    cv::merge(splitResult, returnImg);

    //namedWindow("win", cv::WINDOW_NORMAL);  
    //imshow("win", returnImg);
    //waitKey(0);



    return returnImg;
}
//#endif // UTILS_GENERAL_TIM_H
