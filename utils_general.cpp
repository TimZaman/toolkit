//utils_general.cpp

#include "utils_general.h"

std::vector<std::string> util::split(std::string data, std::string token) {
    std::vector<std::string> output;
    unsigned long pos = std::string::npos;
    do {
        pos = data.find(token);
        output.push_back(data.substr(0, pos));
        if (std::string::npos != pos)
            data = data.substr(pos + token.size());
    } while (std::string::npos != pos);
    return output;
}


bool util::anySubstringInString(std::vector<std::string> & vector_of_substrings, std::string str) {
    for (int i = 0; i < vector_of_substrings.size(); i++) {
        if (str.find(vector_of_substrings[i]) != std::string::npos) {
            return true;
        }
    }
    return false;
}

std::string util::urlencode(const std::string &s) {
    static const char lookup[]= "0123456789abcdef";
    std::stringstream e;
    for(int i=0, ix=s.length(); i<ix; i++){
        const char& c = s[i];
        if ( (48 <= c && c <= 57) ||//0-9
                (65 <= c && c <= 90) ||//abc...xyz
                (97 <= c && c <= 122) || //ABC...XYZ
                (c=='-' || c=='_' || c=='.' || c=='~')  ) {
            e << c;
        } else {
            e << '%';
            e << lookup[ (c&0xF0)>>4 ];
            e << lookup[ (c&0x0F) ];
        }
    }
    return e.str();
}

double util::calcMedian(std::vector<double> scores) { // Calculates median
    double median = 0;
    size_t size = scores.size();
    if (scores.size() == 0) {
        std::cerr << "Warning! Empty vector in calcMedian!" << std::endl;
        return 0;
    }
    sort(scores.begin(), scores.end());
    if (size  % 2 == 0) {
        median = (scores[size / 2 - 1] + scores[size / 2]) / 2;
    } else {
        median = scores[size / 2];
    }
    return median;
}

double util::calcMeanOfQuarterAndThreeQuarterPercentile(std::vector<double> scores){ 
    //Best function name ever. Takes average of 1/4th and 3/4th percentile (median-style)
    double result=0;

    //assert(scores.size()<3);

    size_t size = scores.size();

    if (scores.size()==0) {
        std::cerr << "Warning! Empty vector in calcMedian!" << std::endl;
        return 0;
    }

    sort(scores.begin(), scores.end());

    int index_quarter = floor(scores.size()*0.25);
    int index_three_quarter = floor(scores.size()*0.75);

    result = (scores[index_quarter]+scores[index_three_quarter])*0.5;

    return result;
}

double util::interpolate(double x, std::vector< std::pair<double, double> > &table) {
    const double INF = std::numeric_limits<double>::max();
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

void util::logASL(std::string str_msg){
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
    asl_log(log_client, NULL, ASL_LEVEL_ERR, str_msg.c_str());
    asl_close(log_client);
  #else
    std::cerr << "logASL(" << str_msg << ")" << std::endl;
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

int util::xfilelength(int fd){ //is static in .h
    struct stat sb;
    if (fstat(fd, &sb) < 0){
        return(-1);
    }
    return(sb.st_size);
}


std::string util::ReplaceAll(std::string str, const std::string& from, const std::string& to) {
    size_t start_pos = 0;
    while((start_pos = str.find(from, start_pos)) != std::string::npos) {
        str.replace(start_pos, from.length(), to);
        start_pos += to.length(); // Handles case where 'to' is a substring of 'from'
    }
    return str;
}

#ifdef UTILS_WITH_BOOST

bool util::isValidURL(std::string strUrl){
    try
    {
        boost::regex re("(ftp|http|https):\/\/.[0-9a-zA-Z:\/.%_-]*");
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

std::vector<std::string> util::regexReplaceInVector(std::vector<std::string> vecNames, std::string strInsert, std::string strRegex){
    std::cout << "regexReplaceInVector([";
    for (int i=0; i<vecNames.size(); i++){ std::cout << vecNames[i] << "], ";}
    std::cout << strInsert << ", " << strRegex << std::endl;    

    boost::regex regexNow(strRegex);
    for (int i=0;i<vecNames.size();i++){
        vecNames[i] = boost::regex_replace(vecNames[i], regexNow, strInsert);
    }
    return vecNames;
}

std::vector<std::string> util::getRegexMatches(std::string str_regex, std::string str_matches){
    std::cout << "getRegexMatches(" << str_regex << "," << str_matches << ")" << std::endl;
    std::vector<std::string> vec_matches;
    boost::regex e(str_regex); 
    boost::match_results<std::string::const_iterator> results;
    if (boost::regex_match(str_matches, results, e)){
        for (int i=1; i<results.size(); i++){
            vec_matches.push_back(results[i]);
        }
    }
    for (int i=0; i<vec_matches.size(); i++){
        std::cout << "vec_matches[" << i << "]=" << vec_matches[i] << std::endl;
    }
    return vec_matches;
}

void util::splitDoubleRegex(std::vector<std::string> & vecMatchesFormat, std::vector<std::string> & vecMatchesFilename){
    std::cout << "splitDoubleRegex([";
    for (int i=0; i<vecMatchesFormat.size(); i++){ std::cout << vecMatchesFormat[i] << ",";}
    std::cout <<"],[";
    for (int i=0; i<vecMatchesFilename.size(); i++){ std::cout << vecMatchesFilename[i] << ",";}
    std::cout << "])" << std::endl;

    boost::regex regex_brace("(\\{.*?\\})(\\{.*?\\})"); //(?<!\\d){(.*?)}
    std::vector<std::string> str_doubles_format(2), str_doubles_name(2);
    for (int i=0; i<vecMatchesFormat.size(); i++){
        //Find the ..}{.. split location
        std::size_t splitpos = vecMatchesFormat[i].find("}{");
        if (splitpos!=std::string::npos){
            std::cout << "Split location found at: " << splitpos << std::endl;
        } else {
            std::cout << "No '{}{}' regex to split found." << std::endl;
            continue;
        }
        splitpos++;
        str_doubles_format[0] = vecMatchesFormat[i].substr(0, splitpos);
        str_doubles_format[1] = vecMatchesFormat[i].substr(splitpos, std::string::npos);


        //correlate the split format to the filename part it corresponds to
        std::string namepart = vecMatchesFilename[i];
        //this is for {digit}{something}
        if (str_doubles_format[0][1]=='%' && str_doubles_format[0][str_doubles_format[0].length()-2]=='d'){
            //This is the digit
            int digitloc=-1;
            for (int id=0; id<namepart.length(); id++){
                if (isdigit(namepart[id])){
                    digitloc=id    ;
                } else {
                    break;
                }
            }
            digitloc++;
            //Split at the point where the digits stop
            str_doubles_name[0] = vecMatchesFilename[i].substr(0, digitloc);
            str_doubles_name[1] = vecMatchesFilename[i].substr(digitloc, std::string::npos);
        }

        //this is for {something}{digit}..


        vecMatchesFilename[i] = str_doubles_name[0];
        vecMatchesFilename.insert( vecMatchesFilename.begin()+i, str_doubles_name[1] );
        vecMatchesFormat[i] = str_doubles_format[0];
        vecMatchesFormat.insert( vecMatchesFormat.begin()+i, str_doubles_format[1] );
        i++; //Since we added another index we need to doubleiterate here
    }

    std::cout << "str_doubles_name[0]=" << str_doubles_name[0] << std::endl;
    std::cout << "str_doubles_name[1]=" << str_doubles_name[1] << std::endl;
    std::cout << "str_doubles_format[0]=" << str_doubles_format[0] << std::endl;
    std::cout << "str_doubles_format[1]=" << str_doubles_format[1] << std::endl;

    std::cout << "end splitDoubleRegex()" << std::endl;
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
    strFormatRx = boost::regex_replace(strFormatEsc, regex_num, "(.*?)"); // @TODO(tzaman) get correct amount of numbers ->THIS CRASHES ON UBUNTU+NLS?
    std::cout << "strFormatRx=" << strFormatRx << std::endl;

    std::vector<std::string> vecMatchesFormat = getRegexMatches(strFormatRx, strFormat);
    std::vector<std::string> vecMatchesFilename = getRegexMatches(strFormatRx, strFilename);

    if (vecMatchesFormat.size() != vecMatchesFilename.size()){
        std::cout << "vecMatches sizes do not correspond, returning." << std::endl;
        return mapFormatAndFile; //No match here.
    }

    //Now here it is still possible to have connected regex string, fe {%06d}{datamatrix}.
    //We will now separate these and give them their own index in the vectors
    splitDoubleRegex(vecMatchesFormat, vecMatchesFilename);

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

std::string util::regex_escape(const std::string& string_to_escape) {
    static const boost::regex re_boostRegexEscape( "[\\^\\.\\$\\|\\(\\)\\[\\]\\*\\+\\?\\/\\\\]" );
    const std::string rep( "\\\\\\1&" );
    std::string result = regex_replace(string_to_escape, re_boostRegexEscape, rep, boost::match_default | boost::format_sed);
    return result;
}

//Extract the correlation names f.e. Original/{%04d}.tif to Original/0014.tif for all outputs in the vector
std::vector<std::string> util::correlateFileVectorFormat(std::vector<std::string> vecFormats, std::string filename, int numAdd, int &numNow, std::vector<std::string> &vecNumFormats){
    std::cout << "correlateFileVectorFormat([";
    for (int i=0; i<vecFormats.size();i++){
        std::cout << vecFormats[i] << ",";
    }
    std::cout << "]," << filename << "," << numAdd << "," << numNow << ")" << std::endl;
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
    
        std::cout << std::string(i->first) << " -> " << std::string(i->second) << std::endl;

        if (i->first.length()<3){
            continue;
        }

        //std::cout << boost::format("%d = %s\n") % i->first % i->second;
        if (i->first[1]=='%' && i->first[i->first.length()-2]=='d'){ //Skip .%*d. files (numberformats)
            number_string = i->second;
            continue;
        }
        //
        vecNumFormats = regexReplaceInVector(vecNumFormats, i->second, escapeRegex(i->first));
    }

    for (int i=0; i<vecNumFormats.size(); i++){
        std::cout << "vecNumFormats[" << i << "]=" << vecNumFormats[i] << std::endl;
    }


    //Now fill in the extracted variables in the vecFILENAMING vector.
    //vecNumFormats = regexReplaceInVector(vecFormats, barcode_now, "\\{datamatrix\\}");

    //string number_string  = getNumberFromFileFormat(vecNumFormats[formatID], filename);
    //std::cout << "number_string=" << number_string << std::endl;

    numNow = std::stoi(number_string) + numAdd;

    if (numAdd != 0) {
        std::string strSize = std::to_string(number_string.size());
        std::string strNumFormat = "%0"+ strSize +"i";  //Makes the format the same size as string number
        number_string = boost::str(boost::format(strNumFormat.c_str()) % (numNow));
    }

    //Now fill in the extracted variables in the vecFILENAMING vector.
    vecFormats = regexReplaceInVector(vecNumFormats, number_string, "\\{(%[0-9]*d)\\}");

    //std::cout << " === SUCCESS correlateFileVectorFormat()" << std::endl;
    return vecFormats;
}

std::string util::changeFileExtension(std::string filename, std::string new_extension) {
    // Note the new extension is supposed to have a period in it: i.e. '.jpg'
    boost::filesystem::path p(filename);
    return p.replace_extension(new_extension).string();
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

#endif // UTILS_WITH_BOOST


#ifdef __EXIF_DATA_H__

// Get an existing tag, or create one if it doesn't exist 
ExifEntry* util::init_tag(ExifData *exif, ExifIfd ifd, ExifTag tag) {
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
ExifEntry* util::create_tag(ExifData *exif, ExifIfd ifd, ExifTag tag, size_t len) {
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

#endif // __EXIF_DATA_H__
