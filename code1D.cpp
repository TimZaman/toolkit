
#include "code1D.h"

using namespace std;
using namespace cv;


std::vector<stripeCode> bc1D::readStripeCode(cv::Mat matImage, double dpi){
	cout << "readStripeCode()" << endl;
	std::vector<stripeCode> vecStripecodes;

	//bool debugstripecode = false;
	bool debugstripecode = false;
	bool useAdaptiveThersholding = true;

	Mat matImageK;
	cvtColor(matImage,matImageK, cv::COLOR_BGR2GRAY);
	cv::Mat matThres;

	// VARIABLES //
	double bar_height_mm_min = 7.2-0.8; //[7.5mm=our NMNH c39] [10.7mm=NMNH cover c39]
	double bar_height_mm_max = 10.7+1.0;

	double bar_ar_min = 11;
	double bar_ar_max = 110;

	int min_characters = 7; //minimum characters in barcode string
	//TODO int max_characters = 20; //maximum characters in barcode string

	double bar_dist_group_mm_max = 15; //Maximum distance between any grouped bar to be part of the bar group

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
		int AT_blocksize = bar_height_px_min*0.15;
		int AT_iseven=AT_blocksize%2;
		AT_blocksize += 1+AT_iseven; //Makes sure the blocksize is an even number
		cout << "AT_blocksize=" << AT_blocksize << endl;
		adaptiveThreshold(matImageK, matThres, 255, ADAPTIVE_THRESH_MEAN_C, THRESH_BINARY_INV, AT_blocksize, 20);
	} else {
		threshold(matImageK, matThres, 127, 255, THRESH_BINARY_INV);
	}

	if (debugstripecode){
		cout << "dpi=" << dpi << endl;
		imwrite("/Users/tzaman/Desktop/bc/matImage.tif", matImage);
		imwrite("/Users/tzaman/Desktop/bc/matThres.tif", matThres);
	}

	vector< vector<Point> > contours;
	vector<Vec4i> hierarchy;
	findContours( matThres, contours, hierarchy, CV_RETR_LIST, CV_CHAIN_APPROX_SIMPLE );
	
	cout << "contours.size()=" << contours.size() << endl;
	
	if (contours.size()==0){
		string strErr = "No contours found.";
		cout << strErr << endl;
		return vecStripecodes;
	}

	//RANSAC vars
	int min_inliers = (min_characters+2)*5*0.8; //+2 (start&stop), *5 (stripes per char), *0.8 (margin)
	double max_px_dist = bar_height_px_max*0.1; //Height deviation for ransac

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

		//cout << i << " rotRect: sz=" << rotRect.size << " cp=" << rotRect.center << " a=" << rotRect.angle << " ar=" << ar << endl;


		Rect rCrop = boundingRect(contours[i]);
		if (debugstripecode){
			imwrite("/Users/tzaman/Desktop/seg/" + std::to_string(i) +  ".tif", matImageK(rCrop));
		}
		stripeCandidates.push_back(rotRect);
	}


//TODO: OUTPUT CANDIDATESif (debugstripecode){
//TODO: OUTPUT CANDIDATES	imwrite("/Users/tzaman/Desktop/bc/_" + barcode + ".tif", matBarcode2D);
//TODO: OUTPUT CANDIDATES	Mat matBarcodeFull = matImage.clone();
//TODO: OUTPUT CANDIDATES	util::rectangle(matBarcodeFull, rotRectBarcode, Scalar(50,255,50), 5);
//TODO: OUTPUT CANDIDATES	for (int j=0; j<stripesVerified.size(); j++){
//TODO: OUTPUT CANDIDATES		circle(matBarcodeFull, stripesVerified[j].center, 5, Scalar(50,50,255), 1, CV_AA,0);
//TODO: OUTPUT CANDIDATES	}
//TODO: OUTPUT CANDIDATES	imwrite("/Users/tzaman/Desktop/bc/_matBarcodeFull.tif", matBarcodeFull);
//TODO: OUTPUT CANDIDATES}


	//cout << "stripeCandidates.size()=" << stripeCandidates.size() << endl;

	if (stripeCandidates.size() < min_inliers){
		string strErr = "Code 39 did not find enough bars to accurately make a code.";
		cout << strErr << endl;
		return vecStripecodes;
	}
	
	std::vector<Point> vecPtRectCenter = util::vecrotrect2vecpt(stripeCandidates);
	std::vector<std::vector<int> > vecGroupIdxs = util::groupPoints(vecPtRectCenter, bar_dist_group_px_max, min_inliers);
	std::vector<std::vector<cv::Point> > vecGroupPts(vecGroupIdxs.size());

	//Relate indexes to points and add to group vector
	for (int i=0; i<vecGroupIdxs.size(); i++){
		vecGroupPts[i].resize(vecGroupIdxs[i].size());
		for (int j=0; j<vecGroupIdxs[i].size(); j++){
			//cout << i << "," << j << endl;
			vecGroupPts[i][j] = vecPtRectCenter[vecGroupIdxs[i][j]];
		}
	}

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

	if (vecGroupPts.size()==0){
		string strErr = "Code 39 failed to ransac bars in a line.";
		cout << strErr << endl;
		return vecStripecodes;
	}

	//Now cycle over the groups
	vector<vector<int> > vecVecInlierIdx;
	vector<Vec4f> vecLines;
	vector<int> vecFromGroup; //Keeps track of which group the vecvecInlierIdx belongs to
	for (int i=0; i<vecGroupPts.size(); i++){
		cpRansac_barcode(vecGroupPts[i], min_inliers, max_px_dist, vecVecInlierIdx, vecLines, matImage);
		vecFromGroup.resize(vecVecInlierIdx.size(), i);
	}

	if (vecLines.size()==0){
		string strErr = "Code 39 failed to ransac bars in a line.";
		cout << strErr << endl;
		return vecStripecodes;
	} else {
		cout << "Code39 ransac succesfull" << endl;
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
			double dim_smallest = min(stripesVerified[j].size.width, stripesVerified[j].size.height);
			double dim_tallest  = max(stripesVerified[j].size.width, stripesVerified[j].size.height);
			bar_heights[j] = dim_tallest;
			bar_widths[j]  = dim_smallest;

			//Rotate the points straight
			Point2f ptRot = util::rotatePoint(stripesVerified[j].center, Point(matImageK.cols, matImageK.rows), angle_rad);
			//cout << ptRot << endl;
			coords_x[j]    = ptRot.x;
		}
		
		double height_median = util::calcMedian(bar_heights);
		double width_mean = util::calcMean(bar_widths);
		cout << "height_median=" << height_median <<" width_mean=" << width_mean << endl;

		//Find the start and end position for reading
		vector<size_t> coords_sorted_index;
		vector<double> coords_x_sorted;
		sort(coords_x, coords_x_sorted, coords_sorted_index);
		//cout << coords_x_sorted[0] << " -> " << coords_x_sorted[coords_x_sorted.size()-1] << endl;

		//Get extrema-stripes
		Point2f pt_stripe_left = stripeCandidates[vecVecInlierIdx[i][coords_sorted_index[0]]].center;
		Point2f pt_stripe_right = stripeCandidates[vecVecInlierIdx[i][coords_sorted_index[coords_sorted_index.size()-1]]].center;
		cout << "pt_stripe_left=" << pt_stripe_left << endl;
		cout << "pt_stripe_right=" << pt_stripe_right << endl;

		Point2f pt_barcode_center = (pt_stripe_left+pt_stripe_right)*0.5;
		cout << "pt_barcode_center=" << pt_barcode_center << endl;

		//Calculate width of the barcode
		double barcode_width = util::pointDist(pt_stripe_left, pt_stripe_right);
		cout << "barcode_width=" << barcode_width << endl;

		//Make the rotated rectangle around the barcode
		RotatedRect rotRectBarcode(pt_barcode_center, Size2f(barcode_width, height_median), angle_deg );

		//Add margin (of a few median widths)
		rotRectBarcode.size += Size2f(width_mean*10, 0);



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
		Mat matBarcode1D;
		cv::reduce(matBarcode2D, matBarcode1D, matBarcode2D.cols < matBarcode2D.rows ? 1 : 0, CV_REDUCE_AVG);

		//Make it twice as wide
		double scale_barcode_for_readout = 2.0;
		if (matBarcode1D.cols > matBarcode1D.rows){ 
			resize(matBarcode1D, matBarcode1D, Size(), scale_barcode_for_readout, 1);
		} else { 
			resize(matBarcode1D, matBarcode1D, Size(), 1, scale_barcode_for_readout);
		}

		//Binarize
		util::autoClipBrighten(matBarcode1D,0.05,0.95);
		threshold(matBarcode1D, matBarcode1D, 127, 255, THRESH_BINARY);

		if(debugstripecode){
			imwrite("/Users/tzaman/Desktop/bc/bc2D_" + std::to_string(i) + ".tif", matBarcode2D);
			imwrite("/Users/tzaman/Desktop/bc/bc1D_" + std::to_string(i) + ".tif", matBarcode1D);
		}
		

		//Walk the code and extract words

		char prevbc=255; //white is start and base color
		char curbc=0;
		int lengthnow=0;
		double width_min = scale_barcode_for_readout*width_mean/5.0; //Take a fair margin, things can get messy with thresholding
		double width_max = scale_barcode_for_readout*width_mean*3.0;

		int width_wide_narrow_division = round(scale_barcode_for_readout*width_mean*1.2); //This is the actual division when something is wide or narrow

		vector<char> words;
		for (int j=0; j<matBarcode1D.rows * matBarcode1D.cols; j++){
			curbc = matBarcode1D.at<char>(j); //Current value (1/0)
			if (curbc!=prevbc){ //1/0 switch
				if ((lengthnow > width_min) && (lengthnow < width_max)) {
					words.push_back(lengthnow > width_wide_narrow_division ? 'w' : 'n');
				} else {
					cout << "word ignored. l=" << lengthnow << endl;
					words.push_back('X');
				}
				lengthnow=0; //Reset length
			}
			prevbc = curbc;
			lengthnow++;
		} 

		std::string bc_string(words.begin(), words.end());

		//We now have a word array in which 'n'=narrow, 'w'=wide, 'X'=very wide/invalid


		cout << "Words(" << words.size() <<"):" << endl;
		for (int j=0; j<words.size(); j++){
			cout << words[j];
		}
		cout << endl;


		//Generate the decodingmap
		map<string, char> decoding = generateDecodingMap();

		//Look for start and stop asterisk
		int bc_start_idx=-1;
		int bc_stop_idx=-1;
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
			//return vecStripecodes;
			continue;
		} else if (bc_stop_idx==-1){
			string strErr = "Code 39 stop asterisk not found.";
			cout << strErr << endl;
			//return vecStripecodes;
			continue;
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
			continue;
		}

		cout << "BARCODE: " << barcode << " raw=" << strRaw << endl;

		stripeCode codeNow;
		codeNow.str = barcode;
		codeNow.rotRect = rotRectBarcode;
		codeNow.barcodeType = 3; //Code39
		codeNow.strRaw = strRaw;
		vecStripecodes.push_back(codeNow);


		if (debugstripecode){
			imwrite("/Users/tzaman/Desktop/bc/_" + barcode + ".tif", matBarcode2D);
			Mat matBarcodeFull = matImage.clone();
			util::rectangle(matBarcodeFull, rotRectBarcode, Scalar(50,255,50), 5);
			for (int j=0; j<stripesVerified.size(); j++){
				circle(matBarcodeFull, stripesVerified[j].center, 5, Scalar(50,50,255), 1, CV_AA,0);
			}
			imwrite("/Users/tzaman/Desktop/bc/_matBarcodeFull.tif", matBarcodeFull);
		}
	}


	cout << "END stripecodeTest()" << endl;
	return vecStripecodes;
}





void bc1D::cpRansac_barcode(std::vector<cv::Point> vecPtsIn, int min_inliers, double max_px_dist, std::vector< std::vector<int> > & vecVecInlierIdx, std::vector<cv::Vec4f> & vecLines, cv::Mat image){
	cout << "cpRansac_barcode().." << endl;

	//This is a custom Ransac function that extracts lotsa lines
	RNG rng;

	bool debugransac = false;

	int minstart; //Amount of randomly chosen points RANSAC will start with
	//int min_inliers; //Minimum amount of inliers for succes
	int numiter; //Amount of iterations
	//double max_px_dist;  //max_px_dist [px] as inlier definition

	minstart=2;
	numiter=10000;
	//max_px_dist=ceil(image.cols*0.003); //low; so will only work for lenses with little distortion
	//min_inliers=8;

	//cout << "max_px_dist=" << max_px_dist << "px" << endl;


	int numpts = vecPtsIn.size();

	//This next array keeps track of which points we have already used
	int incrnumsOK[numpts];
	for (int ii=0;ii<numpts;ii++){
		incrnumsOK[ii]=ii;
	}

	for(int i=0;i<numiter;i++){
		Mat img_debug;
		if (debugransac){
			cout << "#i="<<i<< endl;
			image.copyTo(img_debug);
		}
		vector<Point> vecPtCandidate;
		vector<int> vecPtCandidateIdx;

		int inliers_now=minstart;
		//Make the array [0:1:numpts]
		//This array will keep track of the numbers in vecPtsIn that we use
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
				vecPtCandidate.push_back(vecPtsIn[rn]);
				vecPtCandidateIdx.push_back(rn);
				if (debugransac){
					//cout <<"  pt#"<<ii<<" (x,y)=("<<vecPtsIn[rn].x<<","<<vecPtsIn[rn].y<<")"<<endl;
					circle(img_debug, vecPtsIn[rn], 12, Scalar(255,100,100), 2,CV_AA,0);
					//namedWindow("win", 0);
					//imshow( "win", img_debug );
					//waitKey(10);
				}
			}
		}
		if (tries>5000){break;}
		
		//Fit the line
		Vec4f line;
		Point pt1,pt2;
		fitLine( Mat(vecPtCandidate), line, CV_DIST_L2,0,0.01,0.01);

		//Now get two points on this line
		double dd = sqrt((double)line[0]*line[0] + (double)line[1]*line[1]);
		line[0] /= dd;
		line[1] /= dd;
		double t = (float)(image.cols + image.rows);
		pt1.x = round(line[2] - line[0]*t);
		pt1.y = round(line[3] - line[1]*t);
		pt2.x = round(line[2] + line[0]*t);
		pt2.y = round(line[3] + line[1]*t);


		if (0){//debugransac){
			//Shows the line that passes the angle constrant
			cv::line(img_debug, pt1, pt2, Scalar(0,255,255), 10, CV_AA, 0 );
			namedWindow("win", WINDOW_NORMAL);
			Mat img_tmp;
			resize(img_debug, img_tmp, Size(), 0.05, 0.5);
			imshow( "win", img_tmp );
			waitKey(200);
			//cout << line[0] << " "<< line[1] << " "<< line[2] << " "<< line[3] << endl;
		}

		//For every point not in maybe_inliers we iterate and add it to the set
		for(int ii=0;ii<numpts;ii++){
			if ((incrnums[ii]!=-1) && (incrnumsOK[ii]!=-1)){ //Dont chose points already in the model
				Point pt0=vecPtsIn[ii];
				//Calculate the distance between this point and the line (model)
				double dist=abs((pt2.x-pt1.x)*(pt1.y-pt0.y)-(pt1.x-pt0.x)*(pt2.y-pt1.y));
				dist=dist/(sqrt( pow(double(pt2.x-pt1.x),2) + pow(double(pt2.y-pt1.y),2) ));
				if (dist <= max_px_dist){
					incrnums[ii]=-1;//OK; Include it in the index
					vecPtCandidate.push_back(vecPtsIn[ii]);
					vecPtCandidateIdx.push_back(ii);

					inliers_now++;
					if (debugransac){
						circle(img_debug, vecPtsIn[ii], 10, Scalar(100,255,100), 3,CV_AA,0);
					}
				}
			} else {
				continue; //Continue if point already included
			}

		}

		if(inliers_now>=min_inliers){ //If we have found a succesful line
			cout << "RANSAC found a line.." << endl;

			for(int ii=0;ii<numpts;ii++){ //For all the points
				if (incrnums[ii]==-1){ //If this number is included in our line now
					incrnumsOK[ii]=-1; //Exclude the points from the total array so it cant be chosen next time
				}
			}

			fitLine( Mat(vecPtCandidate), line, CV_DIST_L2,0,0.01,0.01);

			vecLines.push_back(line);
			vector<Point> vecPtsNow;
			vecVecInlierIdx.push_back(vecPtCandidateIdx);

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



