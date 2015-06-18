//utils_general.cpp

#include "utils_general.h"



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
		cerr << "The URL regexp is invalid!" << endl;
		throw(e);
	}
	return true;
}


double util::interpolate(double x, vector< pair<double, double> > &table) {
	const double INF = 1.e100;
	// Assumes that "table" is sorted by .first
	// Check if x is out of bound
	if (x > table.back().first) return INF;
	if (x < table[0].first) return -INF;
	vector<pair<double, double> >::iterator it, it2;
	// INFINITY is defined in math.h in the glibc implementation
	it = lower_bound(table.begin(), table.end(), make_pair(x, -INF));
	// Corner case
	if (it == table.begin()) return it->second;
	it2 = it;
	--it2;
	return it2->second + (it->second - it2->second)*(x - it2->first)/(it->first - it2->first);
}

void util::makeBezier(double gamma, double contrast, int N_SEG, vector<int> &lutX, vector<int> &lutY){
	/*for(int i=0;i<N_SEG;i++){	
		lutX[i] = i;
		lutY[i] = i;
	}*/

	//cout << "inside makeBezier(" << gamma << "," << contrast << "," << N_SEG << ")" << endl;
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

	//cout << "Calculating bezier.." << endl;
	//cubic_bezier(x1, y1, x2, y2, x3, y3, x4, y4, N_SEG, xarr, yarr);

	vector<pair<double, double> > table;
	for (int i=0; i < N_SEG; i++){
		double t = (double)i / (double)(N_SEG-1);

		double a = pow((1.0 - t), 3.0);
		double b = 3.0 * t * pow((1.0 - t), 2.0);
		double c = 3.0 * pow(t, 2.0) * (1.0 - t);
		double d = pow(t, 3.0);

		double x = a * x1 + b * x2 + c * x3 + d * x4;
		double y = a * y1 + b * y2 + c * y3 + d * y4;

		//std::cout << "i" << i << " (x,y)=(" << x << "," << y << ")" << std::endl;
		//std::cout << x << "," << y << std::endl;

		//xarr[i]=x;
		//yarr[i]=y;
		table.push_back(make_pair(x,y));		
		//pts[i][0] = x;
		//pts[i][1] = y;
		//return pts;
	}



	//cout << "Sucesfully constructed bezier.." << endl;

	//cout << "Constructing bezier interpolation table.." << endl;

	//for(int i=0;i<N_SEG;i++){	
		
		//sort(table.begin(), table.end()); //We have a span so sorting is not neccesary
	//}
	//cout << "Interpolation table filled in.."

	//cout << "Interpolating.." << endl;
	int value;	
	for(int i=0;i<N_SEG;i++){	
		lutX[i] = i;
		//value = interp_1(i, xarr, yarr, N_SEG);
		value = interpolate(double(i),table);
		value = (value < 0 )? 0 : value; //Values below zero should be zero
		value = (value > (N_SEG-1) )? N_SEG-1 : value; //values above N_SEG-1 should be N_SEG-1
		lutY[i] = value;
		
		//cout << i << "," << interpolate(double(i),table) << endl; 
	}
	//cout << "Interpolation done, discrete (x) bezier constructed." << endl;
}


void util::makeBezier(double gamma, double contrast, int N_SEG, int lutX[], int lutY[]){
	/*for(int i=0;i<N_SEG;i++){	
		lutX[i] = i;
		lutY[i] = i;
	}*/

	//cout << "inside makeBezier(" << gamma << "," << contrast << "," << N_SEG << ")" << endl;
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

	//cout << "Calculating bezier.." << endl;
	//cubic_bezier(x1, y1, x2, y2, x3, y3, x4, y4, N_SEG, xarr, yarr);

	vector<pair<double, double> > table;
	for (int i=0; i < N_SEG; i++){
		double t = (double)i / (double)(N_SEG-1);

		double a = pow((1.0 - t), 3.0);
		double b = 3.0 * t * pow((1.0 - t), 2.0);
		double c = 3.0 * pow(t, 2.0) * (1.0 - t);
		double d = pow(t, 3.0);

		double x = a * x1 + b * x2 + c * x3 + d * x4;
		double y = a * y1 + b * y2 + c * y3 + d * y4;

		//std::cout << "i" << i << " (x,y)=(" << x << "," << y << ")" << std::endl;
		//std::cout << x << "," << y << std::endl;

		//xarr[i]=x;
		//yarr[i]=y;
		table.push_back(make_pair(x,y));		
		//pts[i][0] = x;
		//pts[i][1] = y;
		//return pts;
	}



	//cout << "Sucesfully constructed bezier.." << endl;

	//cout << "Constructing bezier interpolation table.." << endl;

	//for(int i=0;i<N_SEG;i++){	
		
		//sort(table.begin(), table.end()); //We have a span so sorting is not neccesary
	//}
	//cout << "Interpolation table filled in.."

	//cout << "Interpolating.." << endl;
	int value;	
	for(int i=0;i<N_SEG;i++){	
		lutX[i] = i;
		//value = interp_1(i, xarr, yarr, N_SEG);
		value = interpolate(double(i),table);
		value = (value < 0 )? 0 : value; //Values below zero should be zero
		value = (value > (N_SEG-1) )? N_SEG-1 : value; //values above N_SEG-1 should be N_SEG-1
		lutY[i] = value;
		//lutY[i] = i; //FIXME
		//cout << i << "," << interpolate(double(i),table) << endl; 
	}
	//cout << "Interpolation done, discrete (x) bezier constructed." << endl;
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
    cerr << "logASL(" << strMsg << ")" << endl;
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
  if (fstat(fd, &sb) < 0)
    return(-1);
  return(sb.st_size);
}


std::string util::fileformatToRegex(std::string fileformat ){
  //Escape chars for regex first..
  //Change {%04d} to \d{4}
  //Change {datamatrix} to .*
  //Change {manual} to .*

  string strRegex;
  
  if (fileformat.size() < 2){
    cout << "Input fileformat too small (" << fileformat.size() << "). Returning." << endl;
    return strRegex;
  }

  //First escape the regex
  strRegex = regex_escape(fileformat);
  cout << "fileformat=" << fileformat << endl;
  cout << "strRegex=" << strRegex << endl;


  boost::regex regex_num("\\{(%[0-9]*d)\\}");
  strRegex = boost::regex_replace(strRegex, regex_num, "\\\\d{4}");//TODO use good amount of digits
  cout << "strRegex=" << strRegex << endl;

  //boost::regex regex_anyBrace("\\{(.*?)\\}");
  boost::regex regex_anyBrace("(?<!\\\\d)\\{(.*?)\\}");//(?<!\\d){(.*?)}
  strRegex = boost::regex_replace(strRegex, regex_anyBrace, ".*");

  cout << "strRegex=" << strRegex << endl;
  
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

vector<string> util::getRegexMatches(std::string strRegex, std::string strMatches){
  cout << "getRegexMatches(" << strRegex << "," << strMatches << ")" << endl;
  vector<string> vecMatches;
  boost::regex e(strRegex); 
  std::string chat_input(strMatches);
  boost::match_results<std::string::const_iterator> results;
  if (boost::regex_match(chat_input, results, e)){
    for (int i=1; i<results.size(); i++){
      vecMatches.push_back(results[i]);
    }
  }
  for (int i=0;i<vecMatches.size();i++){
    cout << "vecMatches[" << i << "]=" << vecMatches[i] << endl;
  }
  return vecMatches;
}



std::map<std::string, std::string> util::relateFormatAndFile(std::string strFormat, std::string strFilename){
  // ex: relateFormatAndFile("/tif/{var}_{\%04d}.tif","/home/bla/tif/test_01234.tif");
  cout << "relateFormatAndFile(" << strFormat << ", " << strFilename << ")" << endl;
  std::map<std::string, std::string> mapFormatAndFile;

  string strFormatEsc, strFilenameEsc;

  //Escape the strings
  strFormatEsc = escapeRegex(strFormat);

  cout << "Escaped = " << strFormatEsc << ", " << strFilename << "" << endl;


  //    \/tif\/\{var\}_\{%04d\}\.tif
  boost::regex regex_num(".\\{(.*?)\\}");
  string strFormatRx;
  strFormatEsc = ".*"+strFormatEsc; //Match any prefix (like hotfolder dir etc)
  strFormatRx = boost::regex_replace(strFormatEsc, regex_num, "(.*?)"); //TODO get correct amount of numbers    ->THIS CRASHES ON UBUNTU+NLS?
  cout << "strFormatRx=" << strFormatRx << endl;


    vector<string> vecMatchesFormat = getRegexMatches(strFormatRx, strFormat);
    vector<string> vecMatchesFilename = getRegexMatches(strFormatRx, strFilename);

    if (vecMatchesFormat.size() != vecMatchesFilename.size()){
      cout << "vecMatches sizes do not correspond, returning." << endl;
      return mapFormatAndFile; //No match here.
    }

  for (int i=0;i<vecMatchesFormat.size();i++){
    cout << vecMatchesFormat[i] << "=" << vecMatchesFilename[i] << endl;
    mapFormatAndFile[vecMatchesFormat[i]] = vecMatchesFilename[i];
  }

  //cout << " === END relateFormatAndFile()" << endl;
  return mapFormatAndFile;
}



std::vector<std::string> util::correlateFileVectorFormat(std::vector<std::string> vecFormats, std::string filename, int numAdd, int &numNow, std::vector<std::string> &vecNumFormats){
	cout << "correlateFileVectorFormat(.., filename=" << filename << ")" << endl;
	//vecNumFormats is everything in the format replaced except the actual number;
	//numNow is the current number, numAdd is the number to be added.

	std::map<std::string, std::string> mapFormatToFile;
	for (int i=0; i< vecFormats.size();i++){
		cout << vecFormats[i] << endl;
		mapFormatToFile = relateFormatAndFile(vecFormats[i], filename);
		cout << "map size=" << mapFormatToFile.size() << endl;
		//foreach(auto i, mapFormatToFile){
		//  cout << boost::format("%d = %s\n") % i.first % i.second;
		//}
		if (mapFormatToFile.size()>0){
			break; //We found it, only 1 match is possible because we chose 1 filename obviously.
		}
	}

	if (mapFormatToFile.size()==0){
		cout << "mapFormatToFile.size()==0 !" << endl;
		return vecFormats;
	}

	vecNumFormats = vecFormats;
	string number_string;
	/*
	foreach(auto i, mapFormatToFile){
		cout << boost::format("%d = %s\n") % i.first % i.second;
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
	
		cout << boost::format("%d = %s\n") % i->first % i->second;
		if (i->first[1]=='%' && i->first[i->first.length()-2]=='d'){ //Skip .%*d. files (numberformats)
			number_string = i->second;
			continue;
		}
		vecNumFormats = regexReplaceInVector(vecNumFormats, i->second, escapeRegex(i->first));
	}

	for (int i=0; i<vecNumFormats.size(); i++){
		cout << "vecNumFormats[" << i << "]=" << vecNumFormats[i] << endl;
	}


	//Now fill in the extracted variables in the vecFILENAMING vector.
	//vecNumFormats = regexReplaceInVector(vecFormats, barcode_now, "\\{datamatrix\\}");

	//string number_string  = getNumberFromFileFormat(vecNumFormats[formatID], filename);
	cout << "number_string=" << number_string << endl;

	numNow = boost::lexical_cast<int>( number_string ) +numAdd;

	if (numAdd!=0){
		string strSize = boost::lexical_cast<string>(number_string.size());
		string strNumFormat = "%0"+ strSize +"i";  //Makes the format the same size as string number
		number_string = boost::str(boost::format(strNumFormat.c_str()) % (numNow));
	}

	//Now fill in the extracted variables in the vecFILENAMING vector.
	vecFormats = regexReplaceInVector(vecNumFormats, number_string, "\\{(%[0-9]*d)\\}");

	//cout << " === SUCCESS correlateFileVectorFormat()" << endl;
	return vecFormats;
}


std::vector<std::string> util::folderFilesToVector(string folder){
	vector<string> vecFileNames;
	fs::path directory(folder);
	fs::directory_iterator iter(directory), end;
	for(;iter != end; ++iter){
		//if (iter->path().extension() == ".*"){
		string file = iter->path().filename().string();
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
