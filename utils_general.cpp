//utils_general.cpp

#include "utils_general.h"





double util::calcMean(std::vector<double> scores){ //calculates mean
	if (scores.size()==0) {
		std::cerr << "Warning! Empty vector in calcMedian!" << std::endl;
		return 0;
	}
	double sum=0;
	for (int i=0; i < scores.size(); i++){
		sum += scores[i];
	}
	return sum/double(scores.size());
}


double util::calcMedian(std::vector<double> scores){ //calculates median
  double median=0;

  size_t size = scores.size();

  if (scores.size()==0) {
  	std::cerr << "Warning! Empty vector in calcMedian!" << std::endl;
  	return 0;
  }
  
  sort(scores.begin(), scores.end());

  if (size  % 2 == 0){
      median = (scores[size / 2 - 1] + scores[size / 2]) / 2;
  } else {
      median = scores[size / 2];
  }

  return median;
}

bool util::isValidURL(std::string strUrl){
	try
	{
		boost::regex re("(ftp|http|https):\/\/.[0-9a-zA-Z:\/.%]*");
		if (!boost::regex_match(strUrl, re)){
			//throw "Your URL is not formatted correctly!";
			return false;
		} else {
			return true;
		}
	} catch (boost::regex_error& e) {
		std::cerr << "The URL regexp is invalid!" << std::endl;
		throw(e);
	}
	return true;
}


double util::interpolate(double x, std::vector< std::pair<double, double> > &table) {
	const double INF = 1.e100;
	// Assumes that "table" is sorted by .first
	// Check if x is out of bound
	if (x > table.back().first) return INF;
	if (x < table[0].first) return -INF;
	std::vector<std::pair<double, double> >::iterator it, it2;
	// INFINITY is defined in math.h in the glibc implementation
	it = lower_bound(table.begin(), table.end(), std::make_pair(x, -INF));
	// Corner case
	if (it == table.begin()) return it->second;
	it2 = it;
	--it2;
	return it2->second + (it->second - it2->second)*(x - it2->first)/(it->first - it2->first);
}

void util::makeBezier(double gamma, double contrast, int N_SEG, std::vector<int> &lutX, std::vector<int> &lutY){
	/*for(int i=0;i<N_SEG;i++){	
		lutX[i] = i;
		lutY[i] = i;
	}*/

	//std::cout << "inside makeBezier(" << gamma << "," << contrast << "," << N_SEG << ")" << std::endl;
	//Gamma    values are in double from [-1 to 1], 0 being neutral
	//Contrast values are in double from [-1 to 1], 0 being neutral
	//Contrast values increase

	//double xarr[N_SEG];
	//double yarr[N_SEG];

	//Set start point
	double x1=0.0;
	double y1=0.0;
	//Set end point
	double x4=N_SEG-1.0; //fixme is this -1?
	double y4=N_SEG-1.0; //fixme is this -1?

	//We put the 2,3 bezier points on the line between (0,1) and (1,0)
	double x2=0.5*N_SEG*(1.0-gamma + contrast);
	double x3=0.5*N_SEG*(1.0-gamma - contrast);


	//Because we want the two 'in-between' points to lie in the same line ((0,1) and (1,0)), we calculate y2, y3:
	double y2=N_SEG-x2;
	double y3=N_SEG-x3;

	//std::cout << "Calculating bezier.." << std::endl;
	//cubic_bezier(x1, y1, x2, y2, x3, y3, x4, y4, N_SEG, xarr, yarr);

	std::vector< std::pair<double, double> > table;
	for (int i=0; i < N_SEG; i++){
		double t = (double)i / (double)(N_SEG-1);

		double a = pow((1.0 - t), 3.0);
		double b = 3.0 * t * pow((1.0 - t), 2.0);
		double c = 3.0 * pow(t, 2.0) * (1.0 - t);
		double d = pow(t, 3.0);

		double x = a * x1 + b * x2 + c * x3 + d * x4;
		double y = a * y1 + b * y2 + c * y3 + d * y4;

		//std::cout << "i" << i << " (x,y)=(" << x << "," << y << ")" << std::std::endl;
		//std::cout << x << "," << y << std::std::endl;

		//xarr[i]=x;
		//yarr[i]=y;
		table.push_back(std::make_pair(x,y));		
		//pts[i][0] = x;
		//pts[i][1] = y;
		//return pts;
	}



	//std::cout << "Sucesfully constructed bezier.." << std::endl;

	//std::cout << "Constructing bezier interpolation table.." << std::endl;

	//for(int i=0;i<N_SEG;i++){	
		
		//sort(table.begin(), table.end()); //We have a span so sorting is not neccesary
	//}
	//std::cout << "Interpolation table filled in.."

	//std::cout << "Interpolating.." << std::endl;
	int value;	
	for(int i=0;i<N_SEG;i++){	
		lutX[i] = i;
		//value = interp_1(i, xarr, yarr, N_SEG);
		value = interpolate(double(i),table);
		value = (value < 0 )? 0 : value; //Values below zero should be zero
		value = (value > (N_SEG-1) )? N_SEG-1 : value; //values above N_SEG-1 should be N_SEG-1
		lutY[i] = value;
		
		//std::cout << i << "," << interpolate(double(i),table) << std::endl; 
	}
	//std::cout << "Interpolation done, discrete (x) bezier constructed." << std::endl;
}


void util::makeBezier(double gamma, double contrast, int N_SEG, int lutX[], int lutY[]){
	/*for(int i=0;i<N_SEG;i++){	
		lutX[i] = i;
		lutY[i] = i;
	}*/

	//std::cout << "inside makeBezier(" << gamma << "," << contrast << "," << N_SEG << ")" << std::endl;
	//Gamma    values are in double from [-1 to 1], 0 being neutral
	//Contrast values are in double from [-1 to 1], 0 being neutral
	//Contrast values increase

	//double xarr[N_SEG];
	//double yarr[N_SEG];

	//Set start point
	double x1=0.0;
	double y1=0.0;
	//Set end point
	double x4=N_SEG-1.0; //fixme is this -1?
	double y4=N_SEG-1.0; //fixme is this -1?

	//We put the 2,3 bezier points on the line between (0,1) and (1,0)
	double x2=0.5*N_SEG*(1.0-gamma + contrast);
	double x3=0.5*N_SEG*(1.0-gamma - contrast);


	//Because we want the two 'in-between' points to lie in the same line ((0,1) and (1,0)), we calculate y2, y3:
	double y2=N_SEG-x2;
	double y3=N_SEG-x3;

	//std::cout << "Calculating bezier.." << std::endl;
	//cubic_bezier(x1, y1, x2, y2, x3, y3, x4, y4, N_SEG, xarr, yarr);

	std::vector< std::pair<double, double> > table;
	for (int i=0; i < N_SEG; i++){
		double t = (double)i / (double)(N_SEG-1);

		double a = pow((1.0 - t), 3.0);
		double b = 3.0 * t * pow((1.0 - t), 2.0);
		double c = 3.0 * pow(t, 2.0) * (1.0 - t);
		double d = pow(t, 3.0);

		double x = a * x1 + b * x2 + c * x3 + d * x4;
		double y = a * y1 + b * y2 + c * y3 + d * y4;

		//std::cout << "i" << i << " (x,y)=(" << x << "," << y << ")" << std::std::endl;
		//std::cout << x << "," << y << std::std::endl;

		//xarr[i]=x;
		//yarr[i]=y;
		table.push_back(std::make_pair(x,y));		
		//pts[i][0] = x;
		//pts[i][1] = y;
		//return pts;
	}



	//std::cout << "Sucesfully constructed bezier.." << std::endl;

	//std::cout << "Constructing bezier interpolation table.." << std::endl;

	//for(int i=0;i<N_SEG;i++){	
		
		//sort(table.begin(), table.end()); //We have a span so sorting is not neccesary
	//}
	//std::cout << "Interpolation table filled in.."

	//std::cout << "Interpolating.." << std::endl;
	int value;	
	for(int i=0;i<N_SEG;i++){	
		lutX[i] = i;
		//value = interp_1(i, xarr, yarr, N_SEG);
		value = interpolate(double(i),table);
		value = (value < 0 )? 0 : value; //Values below zero should be zero
		value = (value > (N_SEG-1) )? N_SEG-1 : value; //values above N_SEG-1 should be N_SEG-1
		lutY[i] = value;
		//lutY[i] = i; //FIXME
		//std::cout << i << "," << interpolate(double(i),table) << std::endl; 
	}
	//std::cout << "Interpolation done, discrete (x) bezier constructed." << std::endl;
}
















void util::logASL(std::string strMsg){
  #ifdef __APPLE__
    aslclient log_client;
    // ASL_OPT_STDERR     adds stderr as an output file descriptor
    // ASL_OPT_NO_DELAY   connects to the server immediately
    // ASL_OPT_NO_REMOTE  disables remote-control filter adjustment
    log_client = asl_open("Pixel", "The Pixel Facility", ASL_OPT_STDERR);

    //ASL_LEVEL_INFO -> does not get logged?!
    //ASL_LEVEL_NOTICE
    //ASL_LEVEL_ERR -> gets logged, gives no beeps
    //ASL_LEVEL_EMERG -> gets logged, gives 2 big beeps
    //ASL_LEVEL_DEBUG
    asl_log(log_client, NULL, ASL_LEVEL_ERR, strMsg.c_str());
    asl_close(log_client);
  #else
    std::cerr << "logASL(" << strMsg << ")" << std::endl;
  #endif
}




std::string util::escapeRegex(std::string str) {
  str = ReplaceAll(str, "/", "\\/");
  str = ReplaceAll(str, ".", "\\.");
  str = ReplaceAll(str, "{", "\\{");
  str = ReplaceAll(str, "}", "\\}");
  return str;
}

// NLS
// Original/{%04d}.tif
// Bla/Original/1234.tif
// Original\/.*\.tif
// (?<=Original\/)\d{4}(?=.tif)     //obtains only the digits! :)
// SAN
// {datamatrix}/{datamatrix}_{%04d}.tif
// 1823-1234/1823-1234_4321.tif
// .*\/.*_\d{4}\.tif


std::string util::regex_escape(const std::string& string_to_escape) {
	static const boost::regex re_boostRegexEscape( "[\\^\\.\\$\\|\\(\\)\\[\\]\\*\\+\\?\\/\\\\]" );
	const std::string rep( "\\\\\\1&" );
	std::string result = regex_replace(string_to_escape, re_boostRegexEscape, rep, boost::match_default | boost::format_sed);
	return result;
}



int util::xfilelength(int fd){ //is static in .h
	struct stat sb;
	if (fstat(fd, &sb) < 0){
		return(-1);
	}
	return(sb.st_size);
}


std::string util::fileformatToRegex(std::string fileformat ){
	//Escape chars for regex first..
	//Change {%04d} to \d{4}
	//Change {datamatrix} to .*
	//Change {manual} to .*

	std::string strRegex;

	if (fileformat.size() < 2){
	std::cout << "Input fileformat too small (" << fileformat.size() << "). Returning." << std::endl;
	return strRegex;
	}

	//First escape the regex
	strRegex = regex_escape(fileformat);
	std::cout << "fileformat=" << fileformat << std::endl;
	std::cout << "strRegex=" << strRegex << std::endl;


	boost::regex regex_num("\\{(%[0-9]*d)\\}");
	strRegex = boost::regex_replace(strRegex, regex_num, "\\\\d{4}");//TODO use good amount of digits
	std::cout << "strRegex=" << strRegex << std::endl;

	//boost::regex regex_anyBrace("\\{(.*?)\\}");
	boost::regex regex_anyBrace("(?<!\\\\d)\\{(.*?)\\}");//(?<!\\d){(.*?)}
	strRegex = boost::regex_replace(strRegex, regex_anyBrace, ".*");

	std::cout << "strRegex=" << strRegex << std::endl;

	//Now put a a priory regex in
	//strRegex.insert(0, ".*"); //Accept any prependage


	return strRegex;
}



//std::vector<std::string> regexReplaceInVector(vecFILENAMING, number_now, "\\{(%[0-9]*d)\\}"){
std::vector<std::string> util::regexReplaceInVector(std::vector<std::string> vecNames, std::string strInsert, std::string strRegex){
  boost::regex regexNow(strRegex);
  for (int i=0;i<vecNames.size();i++){
    vecNames[i] = boost::regex_replace(vecNames[i], regexNow, strInsert);
  }
  return vecNames;
}





std::string util::ReplaceAll(std::string str, const std::string& from, const std::string& to) {
    size_t start_pos = 0;
    while((start_pos = str.find(from, start_pos)) != std::string::npos) {
        str.replace(start_pos, from.length(), to);
        start_pos += to.length(); // Handles case where 'to' is a substring of 'from'
    }
    return str;
}

std::vector<std::string> util::getRegexMatches(std::string strRegex, std::string strMatches){
	std::cout << "getRegexMatches(" << strRegex << "," << strMatches << ")" << std::endl;
	std::vector<std::string> vecMatches;
	boost::regex e(strRegex); 
	std::string chat_input(strMatches);
	boost::match_results<std::string::const_iterator> results;
	if (boost::regex_match(chat_input, results, e)){
		for (int i=1; i<results.size(); i++){
			vecMatches.push_back(results[i]);
		}
	}
	for (int i=0; i<vecMatches.size(); i++){
		std::cout << "vecMatches[" << i << "]=" << vecMatches[i] << std::endl;
	}
	return vecMatches;
}



std::map<std::string, std::string> util::relateFormatAndFile(std::string strFormat, std::string strFilename){
  // ex: relateFormatAndFile("/tif/{var}_{\%04d}.tif","/home/bla/tif/test_01234.tif");
  std::cout << "relateFormatAndFile(" << strFormat << ", " << strFilename << ")" << std::endl;
  std::map<std::string, std::string> mapFormatAndFile;

  std::string strFormatEsc, strFilenameEsc;

  //Escape the strings
  strFormatEsc = escapeRegex(strFormat);

  std::cout << "Escaped = " << strFormatEsc << ", " << strFilename << "" << std::endl;


  //    \/tif\/\{var\}_\{%04d\}\.tif
  boost::regex regex_num(".\\{(.*?)\\}");
  std::string strFormatRx;
  strFormatEsc = ".*"+strFormatEsc; //Match any prefix (like hotfolder dir etc)
  strFormatRx = boost::regex_replace(strFormatEsc, regex_num, "(.*?)"); //TODO get correct amount of numbers    ->THIS CRASHES ON UBUNTU+NLS?
  std::cout << "strFormatRx=" << strFormatRx << std::endl;


    std::vector<std::string> vecMatchesFormat = getRegexMatches(strFormatRx, strFormat);
    std::vector<std::string> vecMatchesFilename = getRegexMatches(strFormatRx, strFilename);

    if (vecMatchesFormat.size() != vecMatchesFilename.size()){
      std::cout << "vecMatches sizes do not correspond, returning." << std::endl;
      return mapFormatAndFile; //No match here.
    }

  for (int i=0;i<vecMatchesFormat.size();i++){
    std::cout << vecMatchesFormat[i] << "=" << vecMatchesFilename[i] << std::endl;
    mapFormatAndFile[vecMatchesFormat[i]] = vecMatchesFilename[i];
  }

  //std::cout << " === END relateFormatAndFile()" << std::endl;
  return mapFormatAndFile;
}



std::vector<std::string> util::correlateFileVectorFormat(std::vector<std::string> vecFormats, std::string filename, int numAdd, int &numNow, std::vector<std::string> &vecNumFormats){
	std::cout << "correlateFileVectorFormat(.., filename=" << filename << ")" << std::endl;
	//vecNumFormats is everything in the format replaced except the actual number;
	//numNow is the current number, numAdd is the number to be added.

	std::map<std::string, std::string> mapFormatToFile;
	for (int i=0; i< vecFormats.size();i++){
		std::cout << vecFormats[i] << std::endl;
		mapFormatToFile = relateFormatAndFile(vecFormats[i], filename);
		std::cout << "map size=" << mapFormatToFile.size() << std::endl;
		//foreach(auto i, mapFormatToFile){
		//  std::cout << boost::format("%d = %s\n") % i.first % i.second;
		//}
		if (mapFormatToFile.size()>0){
			break; //We found it, only 1 match is possible because we chose 1 filename obviously.
		}
	}

	if (mapFormatToFile.size()==0){
		std::cout << "mapFormatToFile.size()==0 !" << std::endl;
		return vecFormats;
	}

	vecNumFormats = vecFormats;
	std::string number_string;
	/*
	foreach(auto i, mapFormatToFile){
		std::cout << boost::format("%d = %s\n") % i.first % i.second;
		if (i.first[1]=='%' && i.first[i.first.length()-2]=='d'){ //Skip .%*d. files (numberformats)
			number_string = i.second;
			continue;
		}
		vecNumFormats = regexReplaceInVector(vecNumFormats, i.second, escapeRegex(i.first));
	}
	*/




	typedef std::map<std::string, std::string>::iterator it_type;
	for(it_type i = mapFormatToFile.begin(); i != mapFormatToFile.end(); i++) {
		// i->first = key
		// i->second = value
	
		std::cout << boost::format("%d = %s\n") % i->first % i->second;
		if (i->first[1]=='%' && i->first[i->first.length()-2]=='d'){ //Skip .%*d. files (numberformats)
			number_string = i->second;
			continue;
		}
		vecNumFormats = regexReplaceInVector(vecNumFormats, i->second, escapeRegex(i->first));
	}

	for (int i=0; i<vecNumFormats.size(); i++){
		std::cout << "vecNumFormats[" << i << "]=" << vecNumFormats[i] << std::endl;
	}


	//Now fill in the extracted variables in the vecFILENAMING vector.
	//vecNumFormats = regexReplaceInVector(vecFormats, barcode_now, "\\{datamatrix\\}");

	//string number_string  = getNumberFromFileFormat(vecNumFormats[formatID], filename);
	std::cout << "number_string=" << number_string << std::endl;

	numNow = boost::lexical_cast<int>( number_string ) +numAdd;

	if (numAdd!=0){
		std::string strSize = boost::lexical_cast<std::string>(number_string.size());
		std::string strNumFormat = "%0"+ strSize +"i";  //Makes the format the same size as string number
		number_string = boost::str(boost::format(strNumFormat.c_str()) % (numNow));
	}

	//Now fill in the extracted variables in the vecFILENAMING vector.
	vecFormats = regexReplaceInVector(vecNumFormats, number_string, "\\{(%[0-9]*d)\\}");

	//std::cout << " === SUCCESS correlateFileVectorFormat()" << std::endl;
	return vecFormats;
}


std::vector<std::string> util::folderFilesToVector(std::string folder){
	std::vector<std::string> vecFileNames;
	boost::filesystem::path directory(folder);
	boost::filesystem::directory_iterator iter(directory), end;
	for(;iter != end; ++iter){
		//if (iter->path().extension() == ".*"){
		std::string file = iter->path().filename().string();
		vecFileNames.push_back(file);
	}
	return vecFileNames;
}



// Get an existing tag, or create one if it doesn't exist 
ExifEntry* util::init_tag(ExifData *exif, ExifIfd ifd, ExifTag tag){
 	ExifEntry *entry;
 	// Return an existing tag if one exists 
 	if (!((entry = exif_content_get_entry (exif->ifd[ifd], tag)))) {
 	    // Allocate a new entry
 	    entry = exif_entry_new ();
 	    assert(entry != NULL); // catch an out of memory condition 
 	    entry->tag = tag; // tag must be set before calling exif_content_add_entry 

 	    // Attach the ExifEntry to an IFD 
 	    exif_content_add_entry (exif->ifd[ifd], entry);

 	    // Allocate memory for the entry and fill with default data 
 	    exif_entry_initialize (entry, tag);

 	    // Ownership of the ExifEntry has now been passed to the IFD.
 	    // One must be very careful in accessing a structure after
 	    // unref'ing it; in this case, we know "entry" won't be freed
 	    // because the reference count was bumped when it was added to
 	    // the IFD.
 	    //
 	    exif_entry_unref(entry);
 	}
 	return entry;
}

// Create a brand-new tag with a data field of the given length, in the
// given IFD. This is needed when exif_entry_initialize() isn't able to create
// this type of tag itself, or the default data length it creates isn't the
// correct length.
//
ExifEntry* util::create_tag(ExifData *exif, ExifIfd ifd, ExifTag tag, size_t len){
	//ExifContent *content;
 	void *buf;
 	ExifEntry *entry;

 	// Create a memory allocator to manage this ExifEntry 
 	ExifMem *mem = exif_mem_new_default();
 	assert(mem != NULL); // catch an out of memory condition 

 	// Create a new ExifEntry using our allocator
 	entry = exif_entry_new_mem (mem);
 	assert(entry != NULL);

 	// Allocate memory to use for holding the tag data
 	buf = exif_mem_alloc(mem, len);
 	assert(buf != NULL);

 	// Fill in the entry
 	entry->data = (unsigned char*)buf;
 	entry->size = len;
 	entry->tag = tag;
 	entry->components = len;
 	entry->format = EXIF_FORMAT_UNDEFINED;

 	// Attach the ExifEntry to an IFD */
 	exif_content_add_entry (exif->ifd[ifd], entry);

 	//The ExifMem and ExifEntry are now owned elsewhere
	//exif_mem_free(mem, buf); //breaks
 	exif_mem_unref(mem);
 	exif_entry_unref(entry);

 	return entry;
}
