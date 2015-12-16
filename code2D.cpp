
#include "code2D.h"

using namespace std;
using namespace cv;


using zxing::Ref;
using zxing::ArrayRef;
using zxing::LuminanceSource;
using namespace zxing;
using namespace datamatrix;
using namespace zxing::datamatrix;



std::string bc2D::decode_pure_barcode(cv::Mat matImage){
	cout << "bc2D::decode_pure_barcode()" << endl;

	string bcString="";
	Mat matImageK_orig;
	if(matImage.channels() == 3){
		//vector<Mat> imgPlanes(3);
		//split(matImage, imgPlanes);
		//matImageK = imgPlanes[1];//.clone();
		cvtColor(matImage, matImageK_orig, cv::COLOR_BGR2GRAY);
	} else {
		matImageK_orig=matImage;//.clone();
	}

	//Enhange image
	Mat matImageK=matImageK_orig.clone();

	//Auto enhance darkness and brightness
	util::autoClipBrighten(matImageK, 0.10, 0.90);


	//Extract the width and height of the barcode.
	//Currently we only use 10, 12, 14.

	//Sample the grid, maybe do this by resizing the image down to the size of the barcode itself :)
	Mat matBits;
	int sizes[3]={10, 12, 14};

	float stddevs[3];
	for (int i=0; i<3; i++){
		int w = sizes[i];

		resize(matImageK, matBits, Size(w, w), 0, 0, INTER_AREA); 

		Scalar mean,stddev;
		meanStdDev(matBits, mean, stddev);
		stddevs[i]=stddev[0];
		//cout << "mean=" << mean << " stddev=" << stddev << endl;
	}

	//Find highest one
	float max=0;
	int imax=0;
	for (int i=0;i<3;i++){
		if (stddevs[i]>max){
			imax=i;
			max=stddevs[i];
		}
	}

	//We found the width
	int width = sizes[imax];
	int height = width;


	//Reconstruct the image for sampling
	resize(matImageK, matBits, Size(width, height), 0, 0, INTER_AREA); 


	//Finally, check the barcode's orientation, white patch should be in the top right.
	if (matBits.at<uchar>(0, 0) > 128){ //Its in the top left
		util::rot90(matBits,1); //90CW
	} else if (matBits.at<uchar>(height-1, width-1) > 128){ //Its in the bottom right
		util::rot90(matBits,2); //90CCW
	} else if (matBits.at<uchar>(height-1,0) > 128){ //Its in the bottom left
		util::rot90(matBits,3); //180
	} //else: good place already, leave it.


	if (0){
		imwrite("/Users/tzaman/Desktop/bc/matBitsRot.png",matBits);
		cout << "plz press key" << endl;
		char a;
		cin >> a;
	}



	/*
	bool vals[144]= //this is actually an external 12x12 matrix that says 'test'
	{1,0,1,0,1,0,1,0,1,0,1,0,
	 1,0,1,1,0,0,1,1,0,0,1,1,
	 1,1,0,0,1,0,1,0,1,1,1,0,
	 1,1,1,0,0,0,0,1,0,1,0,1,
	 1,0,1,1,0,0,1,0,0,0,1,0,
	 1,1,0,0,1,0,0,0,0,0,0,1,
	 1,1,1,1,0,1,0,0,0,0,0,0,
	 1,0,0,1,0,1,0,1,0,1,0,1,
	 1,1,1,0,1,1,0,1,1,1,1,0,
	 1,0,1,0,1,0,0,0,0,1,0,1,
	 1,1,1,0,0,1,1,1,0,0,1,0,
	 1,1,1,1,1,1,1,1,1,1,1,1};

	bool vals[100]=
		{1,0,1,0,1,0,1,0,1,0,
		 1,1,1,0,0,0,1,0,1,1,
		 1,0,0,0,1,0,1,0,1,0,
		 1,1,0,0,0,0,1,1,1,1,
		 1,0,1,0,1,0,0,0,0,0,
		 1,0,0,0,1,1,0,1,1,1,
		 1,1,0,0,0,0,0,0,0,0,
		 1,1,0,1,1,0,1,1,0,1,
		 1,0,1,0,0,1,1,0,0,0,
		 1,1,1,1,1,1,1,1,1,1};*/


	vector<int> thresholds;
	thresholds.push_back(100); //Attempt a fixed threshold for value 100 [0:255]
	//Find a dynamic threshold for badly printed targets.
	//The strategy used here is finding a threshold so that the entire
	//'connected edge' (bottom left) is filled (black). So find the minimum black
	//value therein
	int minVal = 255;
	for (int x=0;x<matBits.cols;x++){
		int valNow = matBits.at<uchar>(x,matBits.rows-1);
		if (valNow < minVal){
			minVal = valNow;
		}
	}
	for (int y=0; y<matBits.rows; y++){
		int valNow = matBits.at<uchar>(0,y);
		if (valNow < minVal){
			minVal = valNow;
		}
	}
	//cout << "minVal=" << minVal << endl;
	//Now finally add a few points for a little margin.
	minVal = minVal + 5;
	thresholds.push_back(minVal);

	for (int i=0; i<thresholds.size(); i++){
		Ref<BitMatrix> bits(new BitMatrix(width));
		for (int w=0; w<width; w++){
			for (int h=0; h<height; h++){
				if (matBits.at<uchar>(w, h) < thresholds[i]){
					bits->set(h,w); //toggle black=1
				}
			}
		}

		try{
			datamatrix::Decoder decoder_;
			Ref<DecoderResult> decoderResult(decoder_.decode(bits));

			bcString = decoderResult->getText()->getText();

			cout << "Found barcode:" << bcString << endl;
			//cout << decoderResult->getRawBytes() << endl;

		}catch (const zxing::ChecksumException& e) {  
	        cout << "zxing::ChecksumException: " + string(e.what())  << endl; 
		} catch (const zxing::ReaderException& e) {  
	        cout << "zxing::ReaderException: " + string(e.what())  << endl;
	    } catch (const zxing::IllegalArgumentException& e) {  
	        cout << "zxing::IllegalArgumentException: " + string(e.what())  << endl;
	    } catch (const zxing::Exception& e) {  
	        cout << "zxing::Exception: " + string(e.what())  << endl;
	    } catch (const std::exception& e) {  
	        cout << "std::exception: " + string(e.what())  << endl;
	    } catch (...) { //GOTTA CATCH EM ALL *POKEMON*
	    	//POKEMON!
	    	cout << "Pokemon ZXING Catch!" << endl;
	    }
	}

	cout << " decode_pure_barcode END" << endl;
	return bcString;
}







std::string bc2D::decode_image_barcode(const cv::Mat &matImage, vector<int> vecBarcodeTypes, int numwiggles){
	//This is only used for QR codes atm if i am not mistaken
	string bcString="";
	Mat matImageK_orig;
	if(matImage.channels() == 3){
		//vector<Mat> imgPlanes(3);
		//split(matImage, imgPlanes);
		//matImageK = imgPlanes[1];//.clone();
		cvtColor(matImage, matImageK_orig, cv::COLOR_BGR2GRAY);
	} else {
		matImageK_orig=matImage;//.clone();
	}

	for (int w=0; w<numwiggles; w++){
		Mat matImageK = matImageK_orig.clone();
		if (w==0){
			//Nothing for first wiggle

		} else if (w==1){
			//Second wiggle, crop it to the inside
			int bordersz = matImageK.cols*0.1;
			int bordersz_h = floor(bordersz*0.5);
			Rect r(Point(bordersz_h, bordersz_h), Point(matImageK.cols-bordersz, matImageK.rows-bordersz)); //Crop off 7 pixels on all sides
			r = util::constrainRectInSize(r, matImageK.size());
			matImageK = matImageK(r);

		} else if (w==2){
			//Third wiggle, scale down a bit
			resize(matImageK, matImageK, Size(), 0.9, 0.9, INTER_AREA); 

		} else if (w==3){
			//Fourth wiggle, rotate 10 degrees
			util::rotate(matImageK, 10, matImageK);
			//Then crop the inside a bit because rotating introduces black masking edges
			Rect r(Point(5,5), Point(matImageK.cols-10, matImageK.rows-10)); //Crop off 5 pixels on all sides
			matImageK = matImageK(r);
		} else if (w==4){
			//Add border
			int bordersz = matImageK.cols*0.1;
			copyMakeBorder(matImageK, matImageK, bordersz, bordersz, bordersz, bordersz, BORDER_REPLICATE);

		}

		//Copy the data, making sure we have no loose ends or 'smart' data referencing
		Mat matCopy;
		matImageK.copyTo(matCopy);

		try{
			zxingwrapper * zxwrap = new zxingwrapper(matCopy);
			Ref<zxingwrapper> lumisourceRef(zxwrap);

			//Binarizer * binarizer =new GlobalHistogramBinarizer(lumisourceRef);
			Binarizer * binarizer =new HybridBinarizer(lumisourceRef);
			Ref<Binarizer> binarizerRef(binarizer);

			BinaryBitmap * bitmap = new BinaryBitmap(binarizerRef);
			Ref<BinaryBitmap> bitmapRef(bitmap);
			
			DecodeHints hints = DecodeHints(DecodeHints::TRYHARDER_HINT); //memoryleak?
			Ref<MultiFormatReader> mfReader(new MultiFormatReader());

			//Add the different barcode types from vector
			for (int i=0;i<vecBarcodeTypes.size();i++){
				hints.addFormat(zxing::BarcodeFormat::Value(vecBarcodeTypes[i]));
			}


			Ref<Result> result(mfReader->decode(bitmapRef, hints)); 

			bcString=string(result->getText()->getText());
			cout << bcString << endl;
			/*if (result->getResultPoints()->size()==4){
				bc.UL = Point(result->getResultPoints()[0]->getX(), result->getResultPoints()[0]->getY());
				bc.BL = Point(result->getResultPoints()[1]->getX(), result->getResultPoints()[1]->getY());
				bc.UR = Point(result->getResultPoints()[2]->getX(), result->getResultPoints()[2]->getY());
				bc.BR = Point(result->getResultPoints()[3]->getX(), result->getResultPoints()[3]->getY());
				bc.M = (bc.UL + bc.BL +bc.UR + bc.BR)*0.25;
				bc.barcodeType = result->getBarcodeFormat();
			}*/

			//Dont free news or delete anything here, they are all smart pointers that auto delete
		} catch (const zxing::ReaderException& e) {  
	        //cell_result = "zxing::ReaderException: " + string(e.what());  
	    } catch (const zxing::IllegalArgumentException& e) {  
	        //cell_result = "zxing::IllegalArgumentException: " + string(e.what());  
	        util::logASL("IllegalArgumentException Catch:"+string(e.what()));
	    } catch (const zxing::Exception& e) {  
	        //cell_result = "zxing::Exception: " + string(e.what());  
	    } catch (const std::exception& e) {  
	        //cell_result = "std::exception: " + string(e.what());  
	    } catch (...) { //GOTTA CATCH EM ALL *POKEMON*
	    	//POKEMON!
	    	util::logASL("Pokemon ZXING Catch!");
	    }

	    if (!bcString.empty()){
	    	cout << "found code on wiggle:" << w << endl;
	    	break;
	    }
	}

	return bcString;
}





std::string bc2D::readQR(cv::Mat matImage, double dpi){
	cout << "readQR()" << endl;

	//Mat matImage = imread("/Users/tzaman/Desktop/bc.tif");
	//Mat matImage = imread("/Users/tzaman/Desktop/qr.tif");
	//Mat matImage = imread("/Users/tzaman/Desktop/20150505_Oslo_Run/20150505_153650_CAM0.tiff");
	// Mat matImage = imread("/Users/tzaman/Desktop/20150505_Oslo_Run/20150505_153700_CAM0.tiff");
	//Mat matImage = imread("/Users/tzaman/Desktop/20150505_Oslo_Run/20150505_153533_CAM0.tiff");

	//Mat matImage = imread("/Users/tzaman/Desktop/20150505_Oslo_Run/20150505_153518_CAM0.tiff");
	//Mat matImage = imread("/Users/tzaman/Desktop/20150505_Oslo_Run/20150505_153513_CAM0.tiff");
	//Mat matImage = imread("/Users/tzaman/Desktop/20150505_Oslo_Run/20150505_153623_CAM0.tiff");
	//Mat matImage = imread("/Users/tzaman/Desktop/20150505_Oslo_Run/20150505_153639_CAM0.tiff");
	//Mat matImage = imread("/Users/tzaman/Desktop/20150505_Oslo_Run/20150505_153705_CAM0.tiff");
	Mat matImageK; //grayscale image
	Mat matImageC; //contour image

	cvtColor(matImage, matImageK, CV_BGR2GRAY);  //Convert RGB to BGR

	//cout << "blur_start" << endl;
	//Blur it to supress noise and shit
	cv::GaussianBlur(matImageK, matImageK, cv::Size(0, 0), 1.0);
	//cout << "blur_end" << endl;

	threshold(matImageK, matImageC, 130, 255,THRESH_BINARY); //put it high, the featuers are quite black



	vector< vector<Point> > contours; // Vector for storing contour
	vector<Vec4i> hierarchy;
	findContours( matImageC, contours, hierarchy,CV_RETR_CCOMP, CV_CHAIN_APPROX_SIMPLE ); // Find the contours in the image
	
	
	//cout << "contours#=" << contours.size() << endl;
	//namedWindow( "rectangles" );
	//imshow("rectangles", matImage);
	//waitKey(0);

	//Parameters
	//double dpi = 300;
	
	double squareness_threshold = 0.2; //f.e. 0.1=10% deviation aspect ratio width/height
	
	double barcode_width_inch = 0.64;//QR CODE 0.64 inch

	double bc_sizedef_threshold = 0.4; //f.e. 0.1=10% deviation to probable barcode size (1dim)

	//UNUSED - SHOULD USE double bc_margin_extra_inch = 0.05; //add margin to final crop (in inches)

	//double hist_thres = 0.6; //f.e. 0.75=75% histgram comparison (with band-stop histogram filter)

	double qr_square_ratio_dev = 0.35; //allowed size deviation of the *ratios* of the squares areas, fe 0.2=[-20% to +20%]

	double thres_dist_cluster = 2; //threshold of closeness of found centers of multiple square features




	//Histogram correspondance params
	int histSize = 16; //from 0 to 255
	float range[] = { 0, 256 } ; //the upper boundary is exclusive
	const float* histRange = { range };
	bool uniform = true; bool accumulate = false;
	//Define the band-stop filter (that weights the center, gray, low and the whites and blacks high)
	double idealhist[16]={1.0, 1.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0};

	//Caculate barcode dimention from inches and dpi
	double barcode_width_px = dpi*barcode_width_inch;

	//Calculate size to add to rect in pixels
	//Point bc_ptadd(bc_margin_extra_inch*dpi*0.5, bc_margin_extra_inch*dpi*0.5);
	//Size bc_szadd(bc_margin_extra_inch*dpi, bc_margin_extra_inch*dpi);

	//Calculate thresholds
	//Aspect ratio thresholds
	double ar_th_low  = 1.0 - squareness_threshold;
	double ar_th_high = 1.0 + squareness_threshold;

	//Barcode size thresholds
	double bc_th_low  = barcode_width_px * (1-bc_sizedef_threshold);
	double bc_th_high = barcode_width_px * (1+bc_sizedef_threshold);

	// ********************* //
	//        QR CODE        //
	// ********************* //

	

	double aim_ratio_big_mid   = 2.18;
	double aim_ratio_big_small = 4.75;

	

	//Make a square feature
	vector<Point> contour_square;
	contour_square.push_back(Point(0, 0));
	contour_square.push_back(Point(0, 100));
	contour_square.push_back(Point(100, 100));
	contour_square.push_back(Point(100, 0));

	vector<RotatedRect> rectCandidates;

	bool useDebug=false;
	Mat matGui;
	if (useDebug){
		matGui = matImageK.clone();
	}

	//For the QR code, we want to find 3 squares that overlap on the same place
	for( int i = 0; i< contours.size(); i++ ) {// iterate through each contour. 
		RotatedRect rRect = minAreaRect(contours[i]);

		//First check for squareness
		int w = rRect.size.width;
		int h = rRect.size.height; 
		double ar = double(w)/double(h); //Aspect Ratio
		if (ar<ar_th_low || ar>ar_th_high){ //If not square
			continue;
		}

		//For the QR code, test for the smalles square feature to the bigger one in the 3 markers
		int avg = (w+h)*0.5; //Average dimention of the square (w and h avg)
		if (avg<bc_th_low*0.06 || avg> bc_th_high*0.23){
			continue;
		}

		//cout << contours[i] << endl;


		//We now test if the feature matches a square shape indeed
		double shapematch = matchShapes(contours[i], contour_square, CV_CONTOURS_MATCH_I2, 0);
		//cout << "shapematch=" << shapematch << endl;

		if (shapematch>0.02){ //not a square [maximum correct match i've seen so far is 0.02
			continue;
		}

		if (useDebug){
			Mat image = matImageK.clone();
			Point2f vertices[4];
			rRect.points(vertices);
			for (int i = 0; i < 4; i++){
				line(matGui, vertices[i], vertices[(i+1)%4], Scalar(255,0,0), 1);
			}
			//rectangle(matImage, boundRect, Scalar(255,0,0), 1);
			//drawContours( image, contours, 0, Scalar(255), CV_FILLED, 8, hierarchy ); // Draw the largest contour using previously stored index.
			//imshow("rectangles", image);
			//waitKey(0);
			//imwrite(boost::lexical_cast<string>(i)+"_qr_rects.png", image(rRect.boundingRect()));
		}
		
		rectCandidates.push_back(rRect);

	}

	if (useDebug){
		imwrite("0_qr_gui.png", matGui);
	}

	cout << "Found " << rectCandidates.size () << " rectangular candidates" << endl;


	if (rectCandidates.size()<9){
		cout << "Did not find at least 9 candidates, to make 3 markers that make 1 code." << endl;
		return "";
	}


	// *****
	// Shapes to Markers
	// *****


	//vector<int> dist;
	vector< vector <int> > combinemap(rectCandidates.size());
	vector<int> skipmap;
	for (int i=0; i<rectCandidates.size(); i++){
		//cout << "i=" << i << endl;	


		//Check if we already done this one
		bool skip=false;
		for (int s=0; s<skipmap.size(); s++){
			if (i==skipmap[s]){
				skip=true;
				break;
			}
		}
		if (skip){
			continue;
		}

		//Add itself to itself - NOPE
		// combinemap[i].push_back(i); //Not neccesary, it computes the distance to itself already, so it will always add itself
		// skipmap.push_back(i);


		for (int j=0; j<rectCandidates.size(); j++){
			//cout << "j=" << j << endl;		
			//Compute the distance to each other center
			double d = util::pointDist(rectCandidates[i].center, rectCandidates[j].center);
			

			if (d < thres_dist_cluster){
				//cout << "[" << i << "] to [" << j << "]" << " d=" << d << endl;
				combinemap[i].push_back(j);
				skipmap.push_back(j);
			}
		}
	}



	vector< vector<RotatedRect> > qrMarkers;

	for (int i=0; i<combinemap.size(); i++){
		//cout << "combinemap[" << i << "] candidates=" << combinemap[i].size() << endl;

		//If we find 3 candidates, We see if the area ratios add up.
		if (combinemap[i].size()==3){
			vector<double> ratios(3);
			ratios[0] = rectCandidates[combinemap[i][0]].size.area();
			ratios[1] = rectCandidates[combinemap[i][1]].size.area();
			ratios[2] = rectCandidates[combinemap[i][2]].size.area();

			//Sort the ratios from small to large
			vector<size_t> r_index;
			vector<double> r_sorted;
			sort(ratios, r_sorted, r_index);

			//Compute the ratio with the biggest area and the middle, then the smallest
			double r_m = r_sorted[2]/r_sorted[1]; //Idaelly 2.18
			double r_s = r_sorted[2]/r_sorted[0]; //Ideally 4.75
			//cout << "ratios [middle/big]=" << r_m << " [small/big]=" << r_s << endl;

			//Check if the ratios are okay
			//Big to middle first
			if (r_m < aim_ratio_big_mid*(1.0-qr_square_ratio_dev) || r_m > aim_ratio_big_mid*(1.0+qr_square_ratio_dev)){
				cout << "ratios [middle/big] not good:" << r_m << endl;
				continue;
			}
			//Then check big to small
			if (r_s < aim_ratio_big_small*(1.0-qr_square_ratio_dev) || r_s > aim_ratio_big_small*(1.0+qr_square_ratio_dev)){
				cout << "ratios [middle/small] not good:" << r_s << endl;
				continue;
			}

			//TODO: check that the angles are within bounds! rrect.angle, pretty easy, but beware of >45deg issues

			cout << "Found a QR maker." << endl;
			vector<RotatedRect> rects(3);
			rects[0] = rectCandidates[combinemap[i][r_index[0]]]; //Smallest
			rects[1] = rectCandidates[combinemap[i][r_index[1]]]; //Middle
			rects[2] = rectCandidates[combinemap[i][r_index[2]]]; //Biggest

			//rectangle(matImage, rects[0].boundingRect(), Scalar(0,255,0), 1);
			//rectangle(matImage, rects[1].boundingRect(), Scalar(0,255,0), 1);
			//rectangle(matImage, rects[2].boundingRect(), Scalar(0,255,0), 1);

			qrMarkers.push_back(rects);
		}

	}

	cout << "Found qrMakers :" << qrMarkers.size() << endl;

	if (qrMarkers.size()<3){
		cout << "Did not found enough qrMarkers to make up an QR code." << endl;
		return "";
	}

	// *****
	// Markers to Codes
	// *****

	//Now cluster all QR markers, we want 3 markers per QR code.
	//We can compute the maximum marker-to-marker distance from the given qr size in inch
	double max_m2m_dist = barcode_width_inch*dpi*1.13*1.2; //diagonal distance is 1.13x the width of the entire barcode. 1.2 is then the margin



	skipmap.clear(); //Clear it since we used it before
	vector< vector <int> > codecandidates(qrMarkers.size());
	//Loop over the markers and group the markers
	for (int i=0; i<qrMarkers.size(); i++){

		//Check if we already done this one
		bool skip=false;
		for (int s=0; s<skipmap.size(); s++){
			if (i==skipmap[s]){
				skip=true;
				break;
			}
		}
		if (skip){
			continue;
		}

		Point2f ptI = qrMarkers[i][0].center; //Take the center of this marker

		//Compare with other markers
		for (int j=0; j<qrMarkers.size(); j++){
			Point2f ptJ = qrMarkers[j][0].center; //Take the center of this marker
			//Compute distance and see if its smaller than our threshold
			double d = util::pointDist(ptI, ptJ);
			if (d < max_m2m_dist){
				//cout << "[" << i << "] to [" << j << "]" << " d=" << d << endl;
				//Points are close together, probably part of the same qr code.
				codecandidates[i].push_back(j);
				skipmap.push_back(j);
			}
		}
	}



	vector < vector<RotatedRect> > qrCodes;
	for (int i=0; i<codecandidates.size(); i++){
		//cout << "codecandidates[" << i << "] candidates=" << codecandidates[i].size() << endl;

		//If we find 3 candidates, these markers probably make up a QR code
		if (codecandidates[i].size()==3){
			cout << "Found a new code!" << endl;

			//TODO: somewhere rotate the QR code and sort the 3 feature markers

			vector<RotatedRect> rects(3);
			rects[0] = qrMarkers[codecandidates[i][0]][0]; 
			rects[1] = qrMarkers[codecandidates[i][1]][0]; 
			rects[2] = qrMarkers[codecandidates[i][2]][0]; 
			qrCodes.push_back(rects);

			//rectangle(matImage, rects[0].boundingRect(), Scalar(255,0,0), 1);
			//rectangle(matImage, rects[1].boundingRect(), Scalar(255,0,0), 1);
			//rectangle(matImage, rects[2].boundingRect(), Scalar(255,0,0), 1);

			//Now get the enclosing circle around the markers in this code == DIRTY ==, should crop it out and rotate nicely - TODO
			vector<Point2f> ptsQR(3);
			ptsQR[0]= rects[0].center;
			ptsQR[1]= rects[1].center;
			ptsQR[2]= rects[2].center;
			

			Point2f ptCp;
			float radius;
			cv::minEnclosingCircle(ptsQR, ptCp, radius);
			cout << " QR center point ptCp=" << ptCp << " radius=" << radius << endl;

			double sizemulti=1.2;
			//Now crop it out.
			Rect r(Point(ptCp.x-radius*sizemulti, ptCp.y-radius*sizemulti),Size(radius*2*sizemulti,radius*2*sizemulti));

			r = util::constrainRectInSize(r, matImage.size());

			Mat matBarcode = matImage(r);

			//imwrite("QRcode.png", matBarcode);
			
			vector<int> codeType;
			codeType.push_back(12); //6=dmtx, 12=qr;
			string strBarcode = decode_image_barcode(matBarcode, codeType, 5 );

			if (strBarcode.empty()){
				continue;
			}
			return strBarcode;	
		}
	}

	//imwrite("matImage.png", matImage);

	//exit(-1);
	return "";
}





std::string bc2D::readDMTX(cv::Mat matImage, double dpi, double barcode_width_inch_min, double barcode_width_inch_max, int bin_thres){
	cout << "readDMTX()" << endl;

	bool debugbarcode = false; //FOR PRODUCTION PUT TO FALSE

	//Mat matImage = imread("/Users/tzaman/Desktop/bc.tif");
	//Mat matImage = imread("/Users/tzaman/Desktop/qr.tif");
	//Mat matImage = imread("/Users/tzaman/Desktop/20150505_Oslo_Run/20150505_153650_CAM0.tiff");
	// Mat matImage = imread("/Users/tzaman/Desktop/20150505_Oslo_Run/20150505_153700_CAM0.tiff");
	//Mat matImage = imread("/Users/tzaman/Desktop/20150505_Oslo_Run/20150505_153533_CAM0.tiff");

	//Mat matImage = imread("/Users/tzaman/Desktop/20150505_Oslo_Run/20150505_153518_CAM0.tiff");
	//Mat matImage = imread("/Users/tzaman/Desktop/20150505_Oslo_Run/20150505_153513_CAM0.tiff");
	//Mat matImage = imread("/Users/tzaman/Desktop/20150505_Oslo_Run/20150505_153623_CAM0.tiff");
	//Mat matImage = imread("/Users/tzaman/Desktop/20150505_Oslo_Run/20150505_153639_CAM0.tiff");
	//Mat matImage = imread("/Users/tzaman/Desktop/20150505_Oslo_Run/20150505_153705_CAM0.tiff");
	Mat matImageK; //grayscale image
	Mat matImageKblur;
	Mat matImageC; //contour image

	cvtColor(matImage, matImageK, CV_BGR2GRAY);  //Convert RGB to BGR

	//cout << "blur_start" << endl;
	//Blur it to supress noise and shit
	cv::GaussianBlur(matImageK, matImageKblur, cv::Size(0, 0), 1.0);
	//cout << "blur_end" << endl;


	threshold(matImageKblur, matImageC, bin_thres, 255,THRESH_BINARY);

	vector< vector<Point> > contours; // Vector for storing contour
	vector<Vec4i> hierarchy;
	findContours( matImageC, contours, hierarchy,CV_RETR_CCOMP, CV_CHAIN_APPROX_SIMPLE ); // Find the contours in the image

	if (debugbarcode){
		imwrite("/Users/tzaman/Desktop/bc/_matImageC.tif", matImageC);
		//cout << "PREZZEZ KEY" << endl;
		//char a;
		//cin >> a;
	}


	double squareness_threshold = 0.1; //f.e. 0.1=10% deviation aspect ratio width/height
	
	double barcode_width_inch = 0.6; //DMTX oslo 0.48inch dmtx  (1dim)
	//double bc_sizedef_threshold = 0.3; //f.e. 0.1=10% deviation to probable barcode size (1dim)

	double bc_margin_extra_inch = 0.05; //add margin to final crop (in inches)

	double hist_thres = 0.60; //f.e. 0.75=75% histgram comparison (with band-stop histogram filter)

	//Histogram correspondance params
	int histSize = 16; //from 0 to 255
	float range[] = { 0, 256 } ; //the upper boundary is exclusive
	const float* histRange = { range };
	bool uniform = true; bool accumulate = false;
	//Define the band-stop filter (that weights the center, gray, low and the whites and blacks high)
	double idealhist[16]={1.0, 1.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0};

	//Caculate barcode dimention from inches and dpi
	//double barcode_width_px = dpi*barcode_width_inch;

	//Calculate size to add to rect in pixels
	Point bc_ptadd(bc_margin_extra_inch*dpi*0.5, bc_margin_extra_inch*dpi*0.5);
	Size bc_szadd(bc_margin_extra_inch*dpi, bc_margin_extra_inch*dpi);

	//Calculate thresholds
	//Aspect ratio thresholds
	double ar_th_low  = 1.0 - squareness_threshold;
	double ar_th_high = 1.0 + squareness_threshold;

	//Barcode size thresholds
	double bc_th_low  = dpi*barcode_width_inch_min; //barcode_width_px * (1-bc_sizedef_threshold);
	double bc_th_high = dpi*barcode_width_inch_max; //barcode_width_px * (1+bc_sizedef_threshold);

	for( int i = 0; i< contours.size(); i++ ) {// iterate through each contour. 
		double a=contourArea( contours[i],false);  //  Find the area of contour
		if (a<25){continue;}
		//cout << i << endl;

		
		RotatedRect rRect = util::minAreaSquare(contours[i]);

		//First check for squareness
		int w = rRect.size.width;
		int h = rRect.size.height; 
		double ar = double(w)/double(h); //Aspect Ratio
		if (ar<ar_th_low || ar>ar_th_high){ //If not square
			continue;
		}

		//Distinguish by size
		int avg = (w+h)*0.5; //Average dimension of the square (w and h avg)
		if (avg<bc_th_low || avg> bc_th_high){
			continue;
		}

		//Put the barcodes as upright a we found them
		rRect = util::fixRotatedRect(rRect);

		double rotangle = rRect.angle;



		Rect boundRect = rRect.boundingRect();
		boundRect = util::constrainRectInSize(boundRect, matImageK.size());

		if (debugbarcode){
			Mat image = matImage.clone();
			util::rectangle(image, rRect, Scalar(0,255,0), 1);

			//util::rectangle(image, rSquare, Scalar(0,0,255), 1);

			/*for(int j=0;j<convex_hull.size();j++){
				cout << " " << j << " " << convex_hull[j] << endl;
				int j2 = j+1 == convex_hull.size() ? 0 : j+1; //Workaround to connect last to first
				line(image, convex_hull[j], convex_hull[j2], Scalar(255,0,0), 1);
			}*/

			//rectangle(image, boundRect, Scalar(255,0,0), 1);
			//drawContours( matImage, contours, 0, Scalar(255), CV_FILLED, 8, hierarchy ); // Draw the largest contour using previously stored index.
			//imshow("rectangles", image);
			imwrite("/Users/tzaman/Desktop/bc/" + boost::lexical_cast<string>(i) + "_dmtx_rect.png", image(boundRect));
		}


		

		Mat matHist = matImageK(boundRect).clone();
 		

		util::autoClipBrighten(matHist, 0.05, 0.95);

		


		//Check if histogram is any good, should have two distinct peaks.
		Mat hist32;
		calcHist( &matHist, 1, 0, Mat(), hist32, 1, &histSize, &histRange, uniform, accumulate );

		//Normalize such that the area is 1 (not the peaks!)
		double histarea=0;
		for (int i=0;i<16;i++){
			histarea += hist32.at<float>(i);
			//cout << hist32.at<float>(i) << " ";
		}
		//cout << endl << "histarea=" << histarea << endl;

		//Account for the histogram area, set cum to 1.
		for (int i=0;i<16;i++){
			hist32.at<float>(i) *= (1.0/histarea);
		}

		//cout << hist32 << endl;

		//Compare with the ideal histogram and compute score
		//Compute cross correction
		double cumscore=0;
		for (int i=0;i<16;i++){
			cumscore += idealhist[i]*hist32.at<float>(i);
		}

		//cout << "cumscore=" << cumscore << endl;

		if (debugbarcode){
			imwrite("/Users/tzaman/Desktop/bc/" + boost::lexical_cast<string>(i)+"_dmtx_rect_hist_c" + std::to_string((int)round(cumscore*100)) +".png", matHist);
		}

		//A cumscore value of 1 is the maximum that can be attained.
		if (cumscore < hist_thres){
			cout << "histogram cumscore too low (" << cumscore << " / " << hist_thres << "). rejecting candidate " << i <<"." << endl;
			//continue;
		}

		//See if we can read it purely..

		Mat matPureCrop = util::crop(matImage,rRect);
		//imwrite("dmtx_purecrop.png", matPureCrop);
		std::string strBarcodePure = decode_pure_barcode(matPureCrop);

		if (!strBarcodePure.empty()){
			cout << "Wow! Found a pure datamatrix (" << strBarcodePure << ")." << endl;
			return strBarcodePure;
		}




		//Expand the size of the barcode a bit
		Rect bRectPlus = boundRect-bc_ptadd+bc_szadd;
		
		//Finally, we add a margin we crop off later, because the rotation introduces black edge markings.
		int masksz =bRectPlus.width*0.5;
		Point ptadd_mask(masksz*0.5,masksz*0.5);
		Size  szadd_mask(masksz, masksz);
		Rect bRectPlusMargin = bRectPlus - ptadd_mask + szadd_mask;

		bRectPlusMargin = util::constrainRectInSize(bRectPlusMargin, matImage.size());
		Mat matBarcode = matImage(bRectPlusMargin).clone(); //Take the RGB image, as its not blurred etc




		//Rotate it
		util::rotate(matBarcode, rotangle, matBarcode);

		//Crop off the margin we introduced
		Point ptrm_mask = ptadd_mask;
		Size  sz_orig=bRectPlus.size();//retain original size
		Rect  rMask(ptrm_mask, sz_orig);
 
		rMask = util::constrainRectInSize(rMask, matBarcode.size());
		matBarcode = matBarcode(rMask); 

		if (debugbarcode){
			imwrite("/Users/tzaman/Desktop/bc/" + boost::lexical_cast<string>(i)+"_dmtx.png", matBarcode);
		}

		vector<int> codeType;
		codeType.push_back(6); //6=dmtx, 12=qr;

		//cout << "Attempting to read the barcode.." << endl;
		std::string strBarcode = decode_image_barcode(matBarcode, codeType, 5);


		//imwrite(boost::lexical_cast<string>(i)+"_dmtx.png", matBarcode);

		if (strBarcode.empty()){
			continue;
		} else {
			cout << "Found barcode:" << strBarcode << endl;
			return strBarcode;
		}
	}


	cout << "readDMTX() end" << endl;
	return "";
}


