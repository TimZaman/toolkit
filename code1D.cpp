
#include "code1D.h"

using namespace std;
using namespace cv;



stripeCode bc1D::decode_stripes_zxing(cv::Mat matImage){
	//The input matrix needs to be a 1 dimensional array, including the entire barcode with whatever scale, 
	//  but the least amount of clutter to its sides.
	cout << "bc1D::decode_stripes_zxing(..)" << endl;
	//imwrite("/Users/tzaman/Desktop/bc/decode_stripes_zxing.tif", matImage);
	
	stripeCode codeNow;

	//Make sure it's 1D
	if (matImage.cols!=1 && matImage.rows!=1){
		cout << "Warning, the matrix is multidimensional. Need at least 1 dim with size 1." << endl;
		return codeNow;
	}

	string bcString="";
	Mat matImageK_orig;
	if(matImage.channels() == 3){
		cvtColor(matImage, matImageK_orig, cv::COLOR_BGR2GRAY);
	} else {
		matImageK_orig=matImage;
	}

	//Enhange image
	Mat matImageK = matImageK_orig.clone();

	//Auto enhance darkness and brightness
	util::autoClipBrighten(matImageK, 0.10, 0.90);

	cv::Mat matBits = matImageK;
	int length = matBits.cols*matBits.rows;

	zxing::Ref<zxing::BitArray> bits(new zxing::BitArray(length));

	for (int i=0; i<length; i++){
		if (matBits.at<uchar>(i)<100){ //threshold val
			bits->set(i); //toggle black=1
			//cout << "1";
		} //else {
			//cout << "0";
		//}
	}
	//cout << endl;

	std::vector< zxing::Ref<zxing::oned::OneDReader> > readers;
	readers.push_back(zxing::Ref<zxing::oned::OneDReader>(new zxing::oned::Code39Reader()));
	readers.push_back(zxing::Ref<zxing::oned::OneDReader>(new zxing::oned::Code128Reader()));

	//@TODO: ONLY MEISE!!!!!!
	//if (contains(this->projectName, "MEISE")){
		readers.push_back(zxing::Ref<zxing::oned::OneDReader>(new zxing::oned::EAN13Reader()));
	//}

	std::vector<int> barcodetypes;
	barcodetypes.push_back(3);
	barcodetypes.push_back(5);

	//@TODO: ONLY MEISE!!!!!!
	//if (contains(this->projectName, "MEISE")){
		barcodetypes[0] = 15; //swap code39 with UPC-A(15)
	//} 

	//Loop over normal reading direction mode and 180deg rotated
	for (int oidx=0; oidx<=1; oidx++){
		if (oidx==1){ //Looking in reverse direction
			bits->reverse(); // reverse the row and continue
		}

		//Look for multiple stripecode types
		for (int i=0; i<readers.size(); i++){
			try {
				zxing::Ref<zxing::Result> result = readers[i]->decodeRow(0, bits);

				//if (result->getResultPoints()->size()==2){
				//	cout << "(" << result->getResultPoints()[0]->getX() << "," << result->getResultPoints()[0]->getY() << ")" << endl;
				//	cout << "(" << result->getResultPoints()[1]->getX() << "," << result->getResultPoints()[1]->getY() << ")" <<  endl;
				//}

				bcString = result->getText()->getText();

				if (oidx==1){ //Reverse direciton: so we reverse the points too
					zxing::ArrayRef< zxing::Ref<zxing::ResultPoint> > points(result->getResultPoints());
					if (points) {
						points[0] = zxing::Ref<zxing::ResultPoint>(new zxing::oned::OneDResultPoint(length - points[0]->getX() - 1,
						                                                 points[0]->getY()));
						points[1] = zxing::Ref<zxing::ResultPoint>(new zxing::oned::OneDResultPoint(length - points[1]->getX() - 1,
					                                                 points[1]->getY()));
					}
				}

			}catch (const zxing::ChecksumException& e) {  
				//cout << "zxing::ChecksumException: " + string(e.what())  << endl; 
			} catch (const zxing::ReaderException& e) {  
				//cout << "zxing::ReaderException: " + string(e.what())  << endl;
			} catch (const zxing::IllegalArgumentException& e) {  
				//cout << "zxing::IllegalArgumentException: " + string(e.what())  << endl;
			} catch (const zxing::Exception& e) {  
				//cout << "zxing::Exception: " + string(e.what())  << endl;
			} catch (const std::exception& e) {  
				//cout << "std::exception: " + string(e.what())  << endl;
			} catch (...) { //GOTTA CATCH EM ALL *POKEMON*
				//POKEMON!
				//cout << "Pokemon ZXING Catch!" << endl;
			}
			if (!bcString.empty()){
				//Found!
				cout << "bc1D::decode_stripes_zxing() found: " << bcString << endl;

				codeNow.str = bcString;
				//codeNow.strRaw = ?
				//codeNow.rotRect = Rect(Point(result->getResultPoints()[0]->getX(), result->getResultPoints()[0]->getY()), 
				//                       Point(result->getResultPoints()[1]->getX(), result->getResultPoints()[1]->getY()));
				codeNow.barcodeType = barcodetypes[i];

				return codeNow;
			}
		}
	}
	return codeNow;
}



stripeCode bc1D::decode_c39_tzaman(cv::Mat matBarcode1D){
	stripeCode codeNow;
	//Walk the code and extract words

	vector<double> intervals_black;
	vector<double> intervals_white;
	//We will now perform smart measurement of the white and black stripes.
	char prev_char=matBarcode1D.at<char>(0);
	int current_length=0;
	for (int i=1; i<matBarcode1D.rows * matBarcode1D.cols; i++){
		char current_char=matBarcode1D.at<char>(i);
		current_length++;
		if (current_char!=prev_char){
			//Switch
			if (prev_char==0){
				intervals_black.push_back(current_length);
			} else {
				intervals_white.push_back(current_length);
			}
			current_length=0; //reset
		}
		prev_char = current_char;
	}

	intervals_white.erase(intervals_white.begin()); //erase first white index (as the barcode itself starts with black, but background with white)

	int width_wide_narrow_division_white = util::calcMeanOfQuarterAndThreeQuarterPercentile(intervals_white);
	int width_wide_narrow_division_black = util::calcMeanOfQuarterAndThreeQuarterPercentile(intervals_black);
	//Account for the fact that character spacing is thin and more thin divisors that wide
	width_wide_narrow_division_white *= 1.3;

	//Account for the fact that at most 2/5 are black wide.
	width_wide_narrow_division_black *= 1.1;

	double width_max = max(width_wide_narrow_division_white, width_wide_narrow_division_black) * 2.5;

	/*
	cout << "intervals_white=" << endl;
	for (int i=0; i<intervals_white.size(); i++){
		cout << intervals_white[i] << ",";
	}
	cout << endl;
	cout << "intervals_black=" << endl;
	for (int i=0; i<intervals_black.size(); i++){
		cout << intervals_black[i] << ",";
	}
	cout << endl;
	cout << "width_wide_narrow_division_white=" << width_wide_narrow_division_white << endl;
	cout << "width_wide_narrow_division_black=" << width_wide_narrow_division_black << endl;
	*/


	vector<char> words;
	char prevbc=255; //white is start and base color
	int lengthnow=0;
	for (int i=0; i<matBarcode1D.rows * matBarcode1D.cols; i++){
		char curbc = matBarcode1D.at<char>(i); //Current value (1/0)
		if (curbc!=prevbc){ //1/0 switch
			//if ((lengthnow > width_min) && (lengthnow < width_max)) {
			if ((lengthnow < width_max)) {
				int division_width_threshold;
				if (prevbc==0){
					division_width_threshold = width_wide_narrow_division_black;
				} else {
					division_width_threshold = width_wide_narrow_division_white;
				}
				words.push_back(lengthnow > division_width_threshold ? 'w' : 'n');
				
			} else {
				//cout << "word ignored. l=" << lengthnow << endl;
				words.push_back('X');
			}
			//cout << lengthnow << " ";
			lengthnow=0; //Reset length
		}
		prevbc = curbc;
		lengthnow++;
	} 
	//cout << endl;

	std::string bc_string(words.begin(), words.end());

	//We now have a word array in which 'n'=narrow, 'w'=wide, 'X'=very wide/invalid

	/*
	cout << "Words(" << words.size() <<"):" << endl;
	for (int j=0; j<words.size(); j++){
		cout << words[j];
	}
	cout << endl;
	*/


	//Generate the decodingmapc
	map<string, char> decoding = generateDecodingMap();

	//Look for start and stop asterisk
	int bc_start_idx = -1;
	int bc_stop_idx = -1;
	int max_asterisk_start_bars = 9;//Maximum bars that can be passed before asterisk is found

	for (int rev=0; rev<=1;rev++){
		for (int j=0; j<max_asterisk_start_bars; j++){
			map<string, char>::const_iterator it;
			string curr = bc_string.substr(j, 9);
			it = decoding.find(curr);
			if(it == decoding.end()){
				//cout << "NOPE." << endl;
			} else if (it->second==C39_SENTINEL) {
				bc_start_idx=j+10;
				//cout << "FOUND:" << it->second << endl;
			}
		}

		
		for (int j=bc_string.length()-9; j>bc_string.length()-max_asterisk_start_bars-9; j--){
			map<string, char>::const_iterator it;
			string curr = bc_string.substr(j, 9);
			it = decoding.find(curr);
			if(it == decoding.end()){
				//cout << "NOPE." << endl;
			} else if (it->second==C39_SENTINEL) {
				bc_stop_idx=j;
				//cout << "FOUND:" << it->second << endl;
			}
		}
		if (bc_start_idx!=-1 && bc_stop_idx!=-1){ 
			break; //Found!
		} else { //Not found!
			//Try in reverse, reverse the string
			
			std::reverse(std::begin(bc_string), std::end(bc_string));
		}
	}


	if (bc_start_idx==-1){
		string strErr = "Code 39 start asterisk not found.";
		cout << strErr << endl;
		return codeNow;
	} else if (bc_stop_idx==-1){
		string strErr = "Code 39 stop asterisk not found.";
		cout << strErr << endl;
		return codeNow;
	}

	cout << "bc_start_idx=" << bc_start_idx << " bc_stop_idx=" << bc_stop_idx << endl;

	string strRaw = bc_string.substr(bc_start_idx, bc_stop_idx-bc_start_idx+9);

	/*if (((bc_stop_idx-bc_start_idx)%10)!=0){
		string strErr = "Code39 bar count not a modulo of 10, len=" + std::to_string(bc_stop_idx-bc_start_idx);
		cout << strErr << endl;
		return vecStripecodes;
	}*/
	
	//Decode the bitch
	string barcode;
	bool decodingAllOk = true;
	for (int j=bc_start_idx; j<bc_stop_idx; j+=10){
		string curr = bc_string.substr(j, 9);
		//cout << j << " curr=" << curr << endl;

		map<string, char>::const_iterator it;
		it = decoding.find(curr);
		if(it == decoding.end()){
			string strErr = "Code 39 unknown decoding sequence: '" + curr + "' from position " + std::to_string(j); 
			cout << strErr << endl;
			//return vecStripecodes;
			decodingAllOk=false;
			break;
		} else {
			barcode +=(it->second);
		}
	}
	if (!decodingAllOk){
		string strErr = "Code 39 found unknown decoding sequence somewhere, ignoring..";
		return codeNow;
	}

	//cout << "BARCODE: " << barcode << " raw=" << strRaw << endl;


	codeNow.str = barcode;
	//codeNow.rotRect = rotRectBarcode;
	codeNow.barcodeType = 3; //Code39
	codeNow.strRaw = strRaw;
	return codeNow;
}




std::vector<stripeCode> bc1D::readStripeCode(cv::Mat matImage, double dpi){ //warning, 'dpi' variable is overwritten later (400dpi)
	cout << "readStripeCode()" << endl;
	std::vector<stripeCode> vecStripecodes;

	bool debugstripecode = false; //@TODO MAKE SURE TO SET ME TO FALSE IN PRODUCTION
	bool useAdaptiveThersholding = true;

	dpi = 400; //this works well for all scales and sizes..

	Mat matImageK;
	cvtColor(matImage,matImageK, cv::COLOR_BGR2GRAY);
	cv::Mat matThres;

	// VARIABLES //
	double bar_height_mm_min = 3.7; //[7.5mm=our NMNH c39] [10.7mm=NMNH cover c39]
	double bar_height_mm_max = 20;

	double bar_ar_min = 4;
	double bar_ar_max = 110;

	int min_characters = 5; //minimum characters in barcode string

	double bar_dist_group_mm_max = 9.0; //Maximum distance between any grouped bar to be part of the bar group

	// COMPUTE //
	double bar_height_px_min = bar_height_mm_min/25.4*dpi;
	double bar_height_px_max = bar_height_mm_max/25.4*dpi;

	double bar_area_px_min = bar_height_px_min*(bar_height_px_min*1.0/bar_ar_max);

	//Dont allow the area to be less than 1px row
	bar_area_px_min = bar_area_px_min < bar_height_px_min ? bar_height_px_min : bar_area_px_min;

	double bar_area_px_max = bar_height_px_max*(bar_height_px_max*1.0/bar_ar_min);

	double bar_dist_group_px_max = bar_dist_group_mm_max/25.4*dpi;

	if (useAdaptiveThersholding){
		//int AT_blocksize = dpi*0.05; 
		int AT_blocksize = bar_height_px_min*0.5;
		int AT_iseven=AT_blocksize%2;
		AT_blocksize += 1+AT_iseven; //Makes sure the blocksize is an even number
		//cout << "AT_blocksize=" << AT_blocksize << endl;
		adaptiveThreshold(matImageK, matThres, 255, ADAPTIVE_THRESH_MEAN_C, THRESH_BINARY, AT_blocksize, 20);
	} else {
		threshold(matImageK, matThres, 127, 255, THRESH_BINARY_INV);
	}

	if (debugstripecode){
		//cout << "dpi=" << dpi << endl;
		imwrite("/Users/tzaman/Desktop/bc/matImage.tif", matImage);
		imwrite("/Users/tzaman/Desktop/bc/matThres.tif", matThres);
	}

	vector< vector<Point> > contours;
	vector<Vec4i> hierarchy;
	findContours( matThres, contours, hierarchy, CV_RETR_LIST, CV_CHAIN_APPROX_SIMPLE );
	
	//cout << "contours.size()=" << contours.size() << endl;
	
	if (contours.size()==0){
		string strErr = "No contours found.";
		cout << strErr << endl;
		return vecStripecodes;
	}

	//RANSAC vars
	int min_inliers = (min_characters+2)*5*0.75; //+2 (start&stop), *5 (stripes per char), *0.x (margin)
	double max_px_dist = (bar_height_px_min+bar_height_px_max)*0.5*0.05; //Maximum distance from RANSAC line to a point

	vector<RotatedRect> stripeCandidates;
	for(int i = 0; i >= 0; i = hierarchy[i][0] ) {
		double cArea = contourArea( contours[i],false); 
		if (cArea < bar_area_px_min*0.5){
			continue;
		}
		if (cArea > bar_area_px_max){
			continue;
		}
		//cout << "[" << i << "]" << " cArea=" << cArea << endl;

		RotatedRect rotRect= minAreaRect(contours[i]);
		double ar = max(double(rotRect.size.width),double(rotRect.size.height)) / min(double(rotRect.size.width),double(rotRect.size.height));

		if (ar < bar_ar_min){
			continue;
		}
		if (ar > bar_ar_max){
			continue;
		}

		double width = std::min(rotRect.size.width, rotRect.size.height);
		double height = std::max(rotRect.size.width, rotRect.size.height);

		//Check the length
		if (height < bar_height_px_min){
			//Stripe too small
			continue;
		}

		if (height > bar_height_px_max ){
			//Stripe too long
			continue;
		}


		//cout << i << " rotRect: sz=" << rotRect.size << " cp=" << rotRect.center << " a=" << rotRect.angle << " ar=" << ar << endl;


		Rect rCrop = boundingRect(contours[i]);

		//Below parameter is dynamic, plz note
		double min_area_fill=0.15;// = 0.25 ;// 0.4 means 40% of the bounding rectangle of the contour needs to be filled
		//The min_area_fill threshold should be dependent on the width in pixels, because there's more noise in thinner ones
		if (width<3){
			min_area_fill = 0.05;
		} else if (width <5){
			min_area_fill = 0.10;
		}


		//Check if the rectangle is actually filled well
		int fullarea = rCrop.area();

		
		if ( (double(cArea)/double(fullarea)) < min_area_fill){
			continue;
		}

		//cout << i << " fullarea=" << fullarea << " carea=" << cArea << endl;

		if (debugstripecode){
			imwrite("/Users/tzaman/Desktop/seg/" + std::to_string(i) +  ".tif", matImageK(rCrop));
		}
		stripeCandidates.push_back(rotRect);
	}


	if (debugstripecode){
		Mat matBarcodeFull = matImage.clone();
		for (int j=0; j<stripeCandidates.size(); j++){
			util::rectangle(matBarcodeFull, stripeCandidates[j], Scalar(255,0,0), 2);
		}
		imwrite("/Users/tzaman/Desktop/bc/_candidates.tif", matBarcodeFull);
	}


	//cout << "stripeCandidates.size()=" << stripeCandidates.size() << endl;

	if (stripeCandidates.size() < min_inliers){
		string strErr = "Code 39 did not find enough bars to accurately make a code.";
		cout << strErr << endl;
		return vecStripecodes;
	}
	
	std::vector<Point> vecPtRectCenter = util::vecrotrect2vecpt(stripeCandidates);
	std::vector<std::vector<int> > vecGroupIdxs = util::groupPoints(vecPtRectCenter, bar_dist_group_px_max, min_inliers);
	//std::vector<std::vector<cv::Point> > vecGroupPts(vecGroupIdxs.size());
	std::vector<std::vector<cv::RotatedRect> > vecGroupRects(vecGroupIdxs.size());

	//Relate indexes to points and add to group vector
	for (int i=0; i<vecGroupIdxs.size(); i++){
		//vecGroupPts[i].resize(vecGroupIdxs[i].size());
		vecGroupRects[i].resize(vecGroupIdxs[i].size());
		for (int j=0; j<vecGroupIdxs[i].size(); j++){
			//cout << i << "," << j << endl;
			//vecGroupPts[i][j] = vecPtRectCenter[vecGroupIdxs[i][j]];
			vecGroupRects[i][j] = stripeCandidates[vecGroupIdxs[i][j]];
		}
	}


	//Draw all groups
	//if(debugstripecode){
	//	for (int i=0; i<vecGroupPts.size(); i++){
	//		Mat matGroup = matImage.clone();
	//		for (int j=0; j<vecGroupPts[i].size(); j++){
	//			circle(matGroup, vecGroupPts[i][j], 5, Scalar(255,0,255), 1, CV_AA,0);				
	//		}
	//		imwrite("/Users/tzaman/Desktop/bc/_group_" + std::to_string(i) + ".tif", matGroup);
	//	}
	//}
	//exit(-1);
		
		


	//cout << "vecGroupPts.size()=" << vecGroupPts.size() << endl;
	//Erase small groups
	//for (int i=vecGroupPts.size()-1; i>=0; i--){
	//	if (vecGroupPts[i].size() < min_inliers){
	//		//Skipping group, too small.
	//		vecGroupIdxs.erase(vecGroupIdxs.begin()+i);
	//		vecGroupPts.erase(vecGroupPts.begin()+i);
	//	}
	//}
	//cout << "vecGroupPts.size()=" << vecGroupPts.size() << endl;

	if (vecGroupIdxs.size()==0){
		string strErr = "Code 39 failed to ransac bars in a line.";
		cout << strErr << endl;
		return vecStripecodes;
	}

	//Now cycle over the groups
	vector<vector<int> > vecVecInlierIdx;
	vector<Vec4f> vecLines;
	vector<int> vecFromGroup; //Keeps track of which group the vecvecInlierIdx belongs to
	for (int i=0; i<vecGroupRects.size(); i++){
		cpRansac_barcode(vecGroupRects[i], min_inliers, max_px_dist, vecVecInlierIdx, vecLines, matImage);
		vecFromGroup.resize(vecVecInlierIdx.size(), i);
	}

	if (vecLines.size()==0){
		string strErr = "Code 39 failed to ransac bars in a line.";
		cout << strErr << endl;
		return vecStripecodes;
	} else {
		//cout << "Code39 ransac succesfull" << endl;
	}

	//for (int i=0; i<vecGroupIdxs.size(); i++){
	//	cout << "Group " << i << " (" << vecGroupIdxs[i].size() << ") : ";
	//	for (int j=0; j<vecGroupIdxs[i].size(); j++){
	//		cout << vecGroupIdxs[i][j] << " ";
	//	}
	//	cout << endl;
	//}


	//Convert back vecVecInlierIdx to original indices
	for (int i=0; i<vecVecInlierIdx.size(); i++){
		//cout << "vecVecInlierIdx[" << i << "] is from group " << vecFromGroup[i] << endl;
		for (int j=0; j<vecVecInlierIdx[i].size(); j++){
			//cout << " " << vecVecInlierIdx[i][j] << " -> " << vecGroupIdxs[vecFromGroup[i]][vecVecInlierIdx[i][j]] << endl;
			vecVecInlierIdx[i][j] = vecGroupIdxs[vecFromGroup[i]][vecVecInlierIdx[i][j]];
		}
	}

	//cout << "Found " << vecLines.size() << " potential stripecodes." << endl;

	//Loop over all found potential barcodes
	for (int i=0; i < vecLines.size(); i++){
		int numpts = vecVecInlierIdx[i].size();
		cout << "Potential barcode #" << i << " with " << numpts << " points." << endl;

		
		//double angle=atan2(vecLines[i][1],vecLines[i][0])*180/M_PI; //For some reason it clips from [-90,90]
		double angle_rad = atan2(vecLines[i][1],vecLines[i][0]); //For some reason it clips from [-90,90]
		double angle_deg = angle_rad*180.0/M_PI;
		//cout << " angle_deg=" << angle_deg << endl;

		vector<double> bar_heights(numpts);
		vector<double> bar_widths(numpts);
		vector<double> coords_x(numpts);
		//Loop over all found and ransac-verified stripes in this barcode
		vector<cv::RotatedRect> stripesVerified(numpts);
		for (int j=0; j < numpts; j++){
			//cout << vecVecInlierIdx[i][j] << endl;
			//cout << "checking out stripecandidate[" << vecVecInlierIdx[i][j] << "] #" << vecVecInlierIdx[i][j] << endl;
			stripesVerified[j] = stripeCandidates[vecVecInlierIdx[i][j]];
			double dim_smallest = min(stripesVerified[j].size.width, stripesVerified[j].size.height); //For rotation invariance
			double dim_tallest  = max(stripesVerified[j].size.width, stripesVerified[j].size.height); //For rotation invariance
			bar_heights[j] = dim_tallest;
			bar_widths[j]  = dim_smallest;

			//Rotate the points straight
			Point2f ptRot = util::rotatePoint(stripesVerified[j].center, Point(matImageK.cols, matImageK.rows), angle_rad);
			//cout << ptRot << endl;
			coords_x[j]    = ptRot.x;
		}
		
		double height_median = util::calcMedian(bar_heights);
		double width_mean = util::calcMean(bar_widths);
		//cout << "height_median=" << height_median <<" width_mean=" << width_mean << endl;

		//Find the start and end position for reading
		vector<size_t> coords_sorted_index;
		vector<double> coords_x_sorted;
		sort(coords_x, coords_x_sorted, coords_sorted_index);
		//cout << coords_x_sorted[0] << " -> " << coords_x_sorted[coords_x_sorted.size()-1] << endl;

		//Get extrema-stripes
		Point2f pt_stripe_left = stripeCandidates[vecVecInlierIdx[i][coords_sorted_index[0]]].center;
		Point2f pt_stripe_right = stripeCandidates[vecVecInlierIdx[i][coords_sorted_index[coords_sorted_index.size()-1]]].center;
		//cout << "pt_stripe_left=" << pt_stripe_left << endl;
		//cout << "pt_stripe_right=" << pt_stripe_right << endl;

		Point2f pt_barcode_center = (pt_stripe_left+pt_stripe_right)*0.5;
		//cout << "pt_barcode_center=" << pt_barcode_center << endl;

		//Calculate width of the barcode
		double barcode_width = util::pointDist(pt_stripe_left, pt_stripe_right);
		//cout << "barcode_width=" << barcode_width << endl;

		//Make the rotated rectangle around the barcode
		RotatedRect rotRectBarcode(pt_barcode_center, Size2f(barcode_width, height_median), angle_deg );

		//Add margin (of a few median widths)
		rotRectBarcode.size += Size2f(width_mean*7, 0);



		//Extract the barcode itself
		cv::RotatedRect rotRectBarcodeThin = rotRectBarcode;
		//Crop off some margin in thickness because we dont want to collapse the entire barcode.
		if (rotRectBarcodeThin.size.width < rotRectBarcodeThin.size.height){
			rotRectBarcodeThin.size.width *= 0.25;
		} else {
			rotRectBarcodeThin.size.height *= 0.25;
		}
		Mat matBarcode2D = util::crop(matImageK, rotRectBarcodeThin);

		//Collapse the barcode
		Mat matBarcode1D, matBarcode1Dmax;
		cv::reduce(matBarcode2D, matBarcode1D, matBarcode2D.cols < matBarcode2D.rows ? 1 : 0, CV_REDUCE_AVG);
		cv::reduce(matBarcode2D, matBarcode1Dmax, matBarcode2D.cols < matBarcode2D.rows ? 1 : 0, CV_REDUCE_MIN);

		//Make it twice as wide
		double scale_barcode_for_readout = 2.0;
		if (matBarcode1D.cols > matBarcode1D.rows){ 
			resize(matBarcode1D, matBarcode1D, Size(), scale_barcode_for_readout, 1);
		} else { 
			resize(matBarcode1D, matBarcode1D, Size(), 1, scale_barcode_for_readout);
		}

		//Binarize
		util::autoClipBrighten(matBarcode1D,0.06,0.94);
		util::autoClipBrighten(matBarcode1Dmax,0.06,0.94);

		
		if(debugstripecode){
			imwrite("/Users/tzaman/Desktop/bc/bc2D_" + std::to_string(i) + ".tif", matBarcode2D);
			imwrite("/Users/tzaman/Desktop/bc/bc1D_" + std::to_string(i) + ".tif", matBarcode1D);
			imwrite("/Users/tzaman/Desktop/bc/bc1Dmax_" + std::to_string(i) + ".tif", matBarcode1Dmax);
		}
		
		stripeCode codeNow;

		//Test three different thresholds to account for different amounts of ink (little ink/a lot of ink/etc)
		int thresholds[]={165, 150, 180, 120, 100, 200, 215};//THIS LIST IS IN ORDER OF BEST PERFORMING TRESHOLDS!
		//int thresholds[]={150};
		for (int j=0; j<7; j++){
			//cout << "j=" << j << endl;
			Mat matBarcode1Dthres, matBarcode1Dmaxthres;
			threshold(matBarcode1D, matBarcode1Dthres, thresholds[j], 255, THRESH_BINARY);
			threshold(matBarcode1Dmax, matBarcode1Dmaxthres, thresholds[j], 255, THRESH_BINARY);

			if(debugstripecode){
				imwrite("/Users/tzaman/Desktop/bc/bc1D_t" + std::to_string(thresholds[j]) + ".tif", matBarcode1Dthres);
				imwrite("/Users/tzaman/Desktop/bc/bc1D_t" + std::to_string(thresholds[j]) + "m.tif", matBarcode1Dmaxthres);
			}
//cout << "1codeNow.str=" << codeNow.str << endl;
			codeNow = decode_stripes_zxing(matBarcode1Dthres);
//cout << "2codeNow.str=" << codeNow.str << endl;
			if (codeNow.str.empty()){
				//My own function did not work, try zxing's
				codeNow = decode_c39_tzaman(matBarcode1Dthres);
//cout << "3codeNow.str=" << codeNow.str << endl;
				if (!codeNow.str.empty()){
					break;
				}

				//try zxing with only taking the max (for a shitty printer)
				codeNow = decode_stripes_zxing(matBarcode1Dmaxthres);
//cout << "4codeNow.str=" << codeNow.str << endl;
				if (!codeNow.str.empty()){
					break;
				}	
			} else {
				break;
			}
		}
//cout << "5codeNow.str=" << codeNow.str << endl;

		if (codeNow.str.empty()){
			continue; //No barcode here, next..
		}
//cout << "6codeNow.str=" << codeNow.str << endl;




		//Add the actual location to the struct
		codeNow.rotRect = rotRectBarcode;

		//Make sure double barcode readout does not happen i.e.: '10007**10008', so split that up
		std::vector<std::string> split_codes;
		split_codes = util::split(codeNow.str, "**");
		if (split_codes.size()>1){
			std::cout << "Found multiple split barcodes (**) from readout:'" << codeNow.str << "'" << std::endl;
		}
		for (int i=0; i<split_codes.size();i++){
			std::cout << "\t Pushing out (sub)barcode #" << i << "=" << split_codes[i] << std::endl;
			stripeCode code_now_sub = codeNow; // Copy its contents
			code_now_sub.str =  split_codes[i]; // And change the string
			vecStripecodes.push_back(code_now_sub); // And push out
		}

		//TODO: removed duplicate barcodes?


		if (debugstripecode){
			imwrite("/Users/tzaman/Desktop/bc/_" + codeNow.str + ".tif", matBarcode2D);
			Mat matBarcodeFull = matImage.clone();
			util::rectangle(matBarcodeFull, rotRectBarcode, Scalar(50,255,50), 5);
			for (int j=0; j<stripesVerified.size(); j++){
				circle(matBarcodeFull, stripesVerified[j].center, 5, Scalar(50,50,255), 1, CV_AA,0);
			}
			imwrite("/Users/tzaman/Desktop/bc/_matBarcodeFull.tif", matBarcodeFull);
		}
	}


	cout << "END bc1D::readStripeCode() - found #" << vecStripecodes.size() << " barcodes." << endl;
	return vecStripecodes;
}



void bc1D::cpRansac_barcode(std::vector<cv::RotatedRect> vecRectsIn/*vecPtsIn*/, int min_inliers, double max_px_dist, std::vector< std::vector<int> > & vecVecInlierIdx, std::vector<cv::Vec4f> & vecLines, cv::Mat image){
	//cout << "cpRansac_barcode().." << endl;

	//This is a custom Ransac function that extracts lotsa lines
	RNG rng;

	const bool debugransac = false;

	const int minstart = 2; //Amount of randomly chosen points RANSAC will start with
	const int numiter = 10000; //Amount of iterations

	const double stripediff_min = 0.88; // 0.88=88%
	const double stripediff_max = 1.12; // 1.12=12%

	int numpts = vecRectsIn.size();

	//This next array keeps track of which points we have already used
	int incrnumsOK[numpts];
	for (int ii=0;ii<numpts;ii++){
		incrnumsOK[ii]=ii;
	}

	for(int i=0; i<numiter; i++){
		
		Mat img_debug;
		if (debugransac){
			cout << "#i="<<i<< endl;
			image.copyTo(img_debug);
		}
		vector<RotatedRect> vecRectsCandidate;
		vector<int> vecRectsCandidateIdx;

		int inliers_now=minstart;
		//Make the array [0:1:numpts]
		//This array will keep track of the numbers in vecRectsIn that we use
		int incrnums[numpts];
		for (int ii=0;ii<numpts;ii++){
			incrnums[ii]=ii;
		}

		//Count how many we have used already
		int numcandidatesLeft=0;
		for (int ii=0;ii<numpts;ii++){
			if(incrnumsOK[ii]!=-1) numcandidatesLeft++;
		}
		if (numcandidatesLeft<min_inliers){
			break;
		}


		//Select 'minstart' amount of points
		int tries=0;
		for(int ii=0;ii<minstart;ii++){
			tries++;
			//cout<<"  ii="<<ii<<endl;
			int rn = floor(rng.uniform(0., 1.)*numpts);
			//cout << "  rn=" << rn << endl;
			if ((incrnums[rn]==-1) || (incrnumsOK[rn]==-1)){ //Make sure the point we chose is unique
				ii--;  //If the point is not unique (already chosen) we select a new one
				if (tries>5000){break;}
			} else {
				incrnums[rn]=-1;
				vecRectsCandidate.push_back(vecRectsIn[rn]);
				vecRectsCandidateIdx.push_back(rn);
				if (debugransac){
					//cout <<"  pt#"<<ii<<" (x,y)=("<<vecPtsIn[rn].x<<","<<vecPtsIn[rn].y<<")"<<endl;
					circle(img_debug, vecRectsIn[rn].center, 12, Scalar(255,100,100), 2,CV_AA,0);
					//namedWindow("win", 0);
					//imshow( "win", img_debug );
					//waitKey(10);
				}
			}
		}
		if (tries>5000){break;}
		
		//We have now selected a few random candidates, compute the stripeheight from those
		int stripeheight_now=0;
		for (int s=0; s<vecRectsCandidate.size(); s++){
			stripeheight_now += max(vecRectsCandidate[s].size.width, vecRectsCandidate[s].size.height); //Sum it up (later divide for average)
		}
		stripeheight_now = stripeheight_now/vecRectsCandidate.size(); //This is the average
		



		//Fit the line
		Vec4f line;
		Point2f pt1,pt2;
		fitLine( cv::Mat(util::vecrotrect2vecpt(vecRectsCandidate)), line, CV_DIST_L2,0,0.01,0.01);

		//Now get two points on this line
		double dd = sqrt((double)line[0]*line[0] + (double)line[1]*line[1]);
		line[0] /= dd;
		line[1] /= dd;
		double t = (float)(image.cols + image.rows);
		pt1.x = round(line[2] - line[0]*t);
		pt1.y = round(line[3] - line[1]*t);
		pt2.x = round(line[2] + line[0]*t);
		pt2.y = round(line[3] + line[1]*t);

		/*
		if (debugransac){
			//Shows the line that passes the angle constrant
			cv::line(img_debug, pt1, pt2, Scalar(0,255,255), 10, CV_AA, 0 );
			namedWindow("win", WINDOW_NORMAL);
			Mat img_tmp;
			resize(img_debug, img_tmp, Size(), 0.05, 0.5);
			imshow( "win", img_tmp );
			waitKey(200);
			//cout << line[0] << " "<< line[1] << " "<< line[2] << " "<< line[3] << endl;
		}*/
		//For every point not in maybe_inliers we iterate and add it to the set
		for(int ii=0;ii<numpts;ii++){
			if ((incrnums[ii]!=-1) && (incrnumsOK[ii]!=-1)){ //Dont chose points already in the model
				Point2f pt0 = vecRectsIn[ii].center;
				//Calculate the distance between this point and the line (model)
				double dist=abs((pt2.x-pt1.x)*(pt1.y-pt0.y)-(pt1.x-pt0.x)*(pt2.y-pt1.y));
				dist=dist/(sqrt( pow(double(pt2.x-pt1.x),2) + pow(double(pt2.y-pt1.y),2) ));
				double longest_side = max(vecRectsIn[ii].size.width, vecRectsIn[ii].size.height);
				double length_diff = longest_side/stripeheight_now;
				if (dist <= max_px_dist && (length_diff > stripediff_min && length_diff < stripediff_max ) ){
					incrnums[ii]=-1;//OK; Include it in the index
					vecRectsCandidate.push_back(vecRectsIn[ii]);
					vecRectsCandidateIdx.push_back(ii);

					inliers_now++;
					if (debugransac){
						circle(img_debug, vecRectsIn[ii].center, 10, Scalar(100,255,100), 3,CV_AA,0);
					}
				}
			} else {
				continue; //Continue if point already included
			}

		}

		if(inliers_now>=min_inliers){ //If we have found a succesful line
			//cout << "RANSAC found a line.." << endl;

			for(int ii=0;ii<numpts;ii++){ //For all the points
				if (incrnums[ii]==-1){ //If this number is included in our line now
					incrnumsOK[ii]=-1; //Exclude the points from the total array so it cant be chosen next time
				}
			}

			fitLine( cv::Mat( util::vecrotrect2vecpt(vecRectsCandidate)), line, CV_DIST_L2, 0, 0.01, 0.01);

			vecLines.push_back(line);
			//vector<Point> vecPtsNow;
			vecVecInlierIdx.push_back(vecRectsCandidateIdx);

			//circle(img_debug, Point(0,y), img_debug.rows/60, Scalar(100,200,100), 2,CV_AA,0);
			if (debugransac){
				double m=10000;
				//cv::line(img_debug, Point(line[2]-m*line[0], line[3]-m*line[1]), Point(line[2]+m*line[0], line[3]+m*line[1]), Scalar(0,0,255), 6, CV_AA, 0 );
				cv::line(img_debug, pt1, pt2, Scalar(0,0,255), 15, CV_AA, 0 );
				//namedWindow("win", WINDOW_NORMAL);
				imwrite("/Users/tzaman/Desktop/bc/img_debug.tif", img_debug); 
				//exit(-1);
				//Mat img_tmp;
				//resize(img_debug, img_tmp, Size(), 0.05, 0.05);
				//imshow( "win", img_tmp );
				//waitKey(5000);
				//break;
			}
		}

		//cout << endl;

	}
}



